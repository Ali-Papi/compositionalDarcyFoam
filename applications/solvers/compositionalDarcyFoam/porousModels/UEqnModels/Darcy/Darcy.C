/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Darcy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace UEqnModels
{
    defineTypeNameAndDebug(Darcy, 0);
    addToRunTimeSelectionTable(UEqnModel, Darcy, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UEqnModels::Darcy::Darcy
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    UEqnModel(dict, interface),
    interface_(interface.modelCast<UEqnModel, relativePhaseInterface>()),
    KrModel_(KrModel::New(dict, interface_)),
    Tauf
    (
        IOobject
        (
            IOobject::groupName("Tauf", interface.name()),
            interface.fluid().mesh().time().timeName(),
            interface.fluid().mesh()
        ),
        interface.fluid().mesh(),
        dimensionedScalar("zero", dimensionSet(-1, 3, 1, 0, 0), Zero)
    ),
    etaf_
    (
        IOobject
        (
            IOobject::groupName("etaf", interface.name()),
            interface.fluid().mesh().time().timeName(),
            interface.fluid().mesh()
        ),
        interface.fluid().mesh(),
        dimensionedScalar("zero", dimless, Zero)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UEqnModels::Darcy::~Darcy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseInterface&
Foam::UEqnModels::Darcy::interface() const
{
    return interface_;
}


void Foam::UEqnModels::Darcy::correct()
{
    const phaseModel& phase = interface_.relative();
    const phaseModel& porous = interface_.reference();

    Tauf = 1./(linearInterpolate(phase.mu())*porous.invKf());
}


void Foam::UEqnModels::Darcy::update()
{}


void Foam::UEqnModels::Darcy::insertTaus
(
    HashPtrTable<surfaceScalarField>& taus
)
{
    const phaseModel& phase = interface_.relative();
    const volScalarField& alpha = phase;
    const phaseSystem& fluid = interface_.fluid();

    volScalarField S(alpha/fluid.alphaVoid());

    surfaceScalarField Spf(fluid.mvConvP()->interpolate(fluid.phip(), S));

    surfaceScalarField alphaTaup
    (
        KrModel_->Krf(Spf)*Tauf
    );

    // Insert Tau
    taus.insert
    (
        IOobject::groupName("Tauf", phase.name()),
        new surfaceScalarField(Tauf)
    );

    // Insert alphaTaup
    taus.insert
    (
        IOobject::groupName("alphaTaup", phase.name()),
        new surfaceScalarField(alphaTaup)
    );

    // Add to alphaTausp (alpha0Taus)
    *taus["alphaTausp"] += alphaTaup;

    // Insert taup
    taus.insert
    (
        IOobject::groupName("taup", phase.name()),
        new surfaceScalarField
        (
            alphaTaup
          / max(fluid.mvConvP()->interpolate(fluid.phip(), alpha), phase.residualAlpha())
        )
    );

    // update expoent eta
    etaf_ = KrModel_->etaf(Spf);

    surfaceScalarField Sgf(fluid.mvConvG()->interpolate(fluid.phig(), S));

    // Insert alphaTaug
    taus.insert
    (
        IOobject::groupName("alphaTaug", phase.name()),
        new surfaceScalarField
        (
            KrModel_->Krf(Sgf)*Tauf
        )
    );
}


void Foam::UEqnModels::Darcy::addTerms
(
    phaseSystem::eqnBlockTable& eqnBlock,
    PtrList<surfaceScalarField>& phi0s
)
{
    const phaseModel& phase = interface_.relative();
    const volScalarField& alpha = phase;
    const phaseSystem& fluid = interface_.fluid();
    const volScalarField& p = fluid.p();

    const word taupName(IOobject::groupName("taup", phase.name()));
    const word alphaTaupName(IOobject::groupName("alphaTaup", phase.name()));

    tmp<fv::convectionScheme<scalar>> mvConvP(fluid.mvConvP());

    // phase flux by pressure
    surfaceScalarField phipi(-fluid.snGradPSf()*fluid.taus(taupName));

    phipi.internalField() *= etaf_.internalField();

    forAll(phase.U()().boundaryField(), patchi)
    {
        if (phase.U()().boundaryField()[patchi].fixesValue())
        {
            phipi.boundaryField()[patchi] == phase.phi()().boundaryField()[patchi];
        }
    }

    if (phase.pure())
    {
        label Ai = phase.blockIndex()[0];
        const surfaceScalarField& alphaTau = fluid.taus(alphaTaupName); 

        surfaceScalarField phiAlpha0P0(alphaTau*fluid.snGradPSf());

        phiAlpha0P0.internalField() *= etaf_.internalField();

        phi0s[Ai] = -phiAlpha0P0;

        fvScalarMatrix divPhiAlpha
        (
            mvConvP->fvmDiv(phipi, alpha)
          + fvc::div(phiAlpha0P0)
        );

        // alpha in alphaEqn [Ai, Ai] (already set)
        eqnBlock[Ai][Ai] += divPhiAlpha;

        // p in alphaEqn [Ai, p]
        fvScalarMatrix pInAi
        (
         - fvm::laplacian(alphaTau, p)
        );

        if (eqnBlock[Ai].set(0))
        {
            eqnBlock[Ai][0] += pInAi;
        }
        else
        {
            eqnBlock[Ai].set(0, pInAi);
        }

        // alpha in pEqn [p, Ai] (already set)
        eqnBlock[0][Ai] += divPhiAlpha;
    }
    else
    {
        forAll(phase.Y(), compi)
        {
            label Ai = phase.blockIndex()[compi];
            const volScalarField& Y = phase.Y()[compi];

            surfaceScalarField alphaYTau(mvConvP->interpolate(fluid.phip(), Y)*fluid.taus(taupName));

            surfaceScalarField phiAlphaY0P0(alphaYTau*fluid.snGradPSf());

            phiAlphaY0P0.internalField() *= etaf_.internalField();

            phi0s[Ai] = -phiAlphaY0P0;

            fvScalarMatrix divPhiAlphaY
            (
                mvConvP->fvmDiv(phipi, Y)
              + fvc::div(phiAlphaY0P0)
            );

            // alpha in alphaYEqn [Ai, Ai] (already set)
            eqnBlock[Ai][Ai] += divPhiAlphaY;

            // p in alphaYEqn [Ai, p]
            fvScalarMatrix pInAi
            (
             - fvm::laplacian(alphaYTau, p)
            );

            if (eqnBlock[Ai].set(0))
            {
                eqnBlock[Ai][0] += pInAi;
            }
            else
            {
                eqnBlock[Ai].set(0, pInAi);
            }

            // alpha in pEqn [p, Ai] (already set)
            eqnBlock[0][Ai] += divPhiAlphaY;
        }
    }
}


// ************************************************************************* //
