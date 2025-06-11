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

#include "pcMultiRockTabulated.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pcModels
{
    defineTypeNameAndDebug(pcMultiRockTabulated, 0);
    addToRunTimeSelectionTable(pcModel, pcMultiRockTabulated, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pcModels::pcMultiRockTabulated::pcMultiRockTabulated
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    pcModel(dict, interface),
    interface_(interface.modelCast<pcModel, relativePhaseInterface>()),
    KrdPcdSf
    (
        IOobject
        (
            IOobject::groupName("KrdPcdSf", interface.name()),
            interface.fluid().mesh().time().timeName(),
            interface.fluid().mesh()
        ),
        interface.fluid().mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), Zero)
    ),
    KrgradYPcf
    (
        IOobject
        (
            IOobject::groupName("KrgradYPcf", interface.name()),
            interface.fluid().mesh().time().timeName(),
            interface.fluid().mesh()
        ),
        interface.fluid().mesh(),
        dimensionedScalar("zero", dimensionSet(1, 0, -2, 0, 0), Zero)
    ),
    scale_(dict.lookupOrDefault("scale", 1.0)),
    porous_(interface.fluid().lookupPhase(dict.lookup("porousPhase")))
{
    if (!porous_.pure() && porous_.porous())
    {
        KrdPcdSTables.setSize(porous_.species().size());
        KrPcSTables.setSize(porous_.species().size());

        forAll(porous_.species(), i)
        {
            const word& speciei = porous_.species()[i];
            const word KrdPcdSName(IOobject::groupName("KrdPcdS", speciei));
            const word KrPcSName(IOobject::groupName("KrPcS", speciei));

            KrdPcdSTables.set
            (
                i,
                new interpolationTable<scalar>(dict.subDict(KrdPcdSName))
            );

            KrPcSTables.set
            (
                i,
                new interpolationTable<scalar>(dict.subDict(KrPcSName))
            );
        }
    }
    else
    {
        FatalErrorInFunction
            << "Selected porousPhase: " << porous_.name()
            << " is not a porous and multicomponent phase."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pcModels::pcMultiRockTabulated::~pcMultiRockTabulated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseInterface&
Foam::pcModels::pcMultiRockTabulated::interface() const
{
    return interface_;
}


void Foam::pcModels::pcMultiRockTabulated::correct()
{
    const phaseModel& phase = interface_.relative();
    const phaseSystem& fluid = interface_.fluid();
    const fvMesh& mesh = interface_.mesh();
    const volScalarField& alpha = phase;
    const volScalarField& alphaVoid = fluid.alphaVoid();

    volScalarField S(alpha/alphaVoid);
    surfaceScalarField Sf(fvc::interpolate(S));

    forAll(KrdPcdSf, facei)
    {
        KrdPcdSf[facei] = 0.;
        KrgradYPcf[facei] = 0.;
    }

    surfaceScalarField::GeometricBoundaryField& KrdPcdSfb = KrdPcdSf.boundaryField();
    surfaceScalarField::GeometricBoundaryField& KrgradYPcfb = KrgradYPcf.boundaryField();
    surfaceScalarField::GeometricBoundaryField& Sfb = Sf.boundaryField();

    forAll(KrdPcdSfb, patchi)
    {
        forAll(KrdPcdSfb[patchi], facei)
        {
            KrdPcdSfb[patchi][facei] = 0.;
            KrgradYPcfb[patchi][facei] = 0.;
        }
    }

    forAll(porous_.Y(), compi)
    {
        const volScalarField& porous = porous_;
        volScalarField Yj(porous_.Y()[compi]/max(porous, 0.001));
        surfaceScalarField snGradYj(fvc::snGrad(Yj)*mesh.magSf());
        surfaceScalarField Yjf(fvc::interpolate(Yj));

        forAll(KrdPcdSf, facei)
        {
            KrdPcdSf[facei] += scale_*Yjf[facei]*KrdPcdSTables[compi](Sf[facei]);
            KrgradYPcf[facei] += scale_*snGradYj[facei]*KrPcSTables[compi](Sf[facei]);
        }

        surfaceScalarField::GeometricBoundaryField& Yjfb = Yjf.boundaryField();
        surfaceScalarField::GeometricBoundaryField& snGradYjb = snGradYj.boundaryField();

        forAll(KrdPcdSfb, patchi)
        {
            if (!KrdPcdSfb[patchi].coupled()) continue;

            forAll(KrdPcdSfb[patchi], facei)
            {
                KrdPcdSfb[patchi][facei] += scale_*Yjfb[patchi][facei]
                  * KrdPcdSTables[compi](Sfb[patchi][facei]);
                KrgradYPcfb[patchi][facei] += scale_*snGradYjb[patchi][facei]
                  * KrPcSTables[compi](Sfb[patchi][facei]);
            }
        }
    }
}


void Foam::pcModels::pcMultiRockTabulated::update()
{}


void Foam::pcModels::pcMultiRockTabulated::addTerms
(
    phaseSystem::eqnBlockTable& eqnBlock,
    PtrList<surfaceScalarField>& phi0s
)
{
    const phaseModel& phase = interface_.relative();
    const volScalarField& alpha = phase;
    const fvMesh& mesh = interface_.mesh();
    const phaseSystem& fluid = interface_.fluid();
    const volScalarField& alphaVoid = fluid.alphaVoid();

    const label& Ai = phase.blockIndex()[0];

    const word TaufName(IOobject::groupName("Tauf", phase.name()));
    const surfaceScalarField& Tauf = fluid.taus(TaufName);

    surfaceScalarField snGradVoid(fvc::snGrad(alphaVoid)*mesh.magSf());
    surfaceScalarField voidf(fvc::interpolate(alphaVoid));
    surfaceScalarField gamma(Tauf*KrdPcdSf/voidf);
    surfaceScalarField phiGradY(Tauf*KrgradYPcf/voidf);

    fvScalarMatrix pcEqn
    (
        fvm::laplacian(gamma, alpha)
      - fvm::div(snGradVoid*gamma/voidf, alpha)
      + fvm::div(phiGradY, alpha, "div(phiGradY,alpha)")
    );

    // alpha in alphaEqn [Ai, Ai] (already set)
    eqnBlock[Ai][Ai] += pcEqn;

    // alpha in pEqn [p, Ai] (already set)
    eqnBlock[0][Ai] += pcEqn;
}


// ************************************************************************* //
