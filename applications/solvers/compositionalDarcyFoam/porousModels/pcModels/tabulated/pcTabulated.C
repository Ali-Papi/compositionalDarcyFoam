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

#include "pcTabulated.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pcModels
{
    defineTypeNameAndDebug(pcTabulated, 0);
    addToRunTimeSelectionTable(pcModel, pcTabulated, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pcModels::pcTabulated::pcTabulated
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
    KrdPcdSTable(dict),
    scale_(dict.lookupOrDefault("scale", 1.0))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pcModels::pcTabulated::~pcTabulated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseInterface&
Foam::pcModels::pcTabulated::interface() const
{
    return interface_;
}


void Foam::pcModels::pcTabulated::correct()
{
    const phaseModel& phase = interface_.relative();
    const phaseSystem& fluid = interface_.fluid();
    const volScalarField& alpha = phase;
    const volScalarField& alphaVoid = fluid.alphaVoid();

    surfaceScalarField Sf(fvc::interpolate(alpha/alphaVoid));

    forAll(KrdPcdSf, facei)
    {
        KrdPcdSf[facei] = scale_*KrdPcdSTable(Sf[facei]);    
    }

    surfaceScalarField::GeometricBoundaryField& KrdPcdSfb = KrdPcdSf.boundaryField();
    surfaceScalarField::GeometricBoundaryField& Sfb = Sf.boundaryField();

    forAll(KrdPcdSfb, patchi)
    {
        if (!KrdPcdSfb[patchi].coupled()) continue;

        forAll(KrdPcdSfb[patchi], facei)
        {
            KrdPcdSfb[patchi][facei] = scale_*KrdPcdSTable(Sfb[patchi][facei]);
        }
    }
}


void Foam::pcModels::pcTabulated::update()
{}


void Foam::pcModels::pcTabulated::addTerms
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

    fvScalarMatrix pcEqn 
    (
        fvm::laplacian(gamma, alpha)
      - fvm::div(snGradVoid*gamma/voidf, alpha)
    );

    // alpha in alphaEqn [Ai, Ai] (already set)
    eqnBlock[Ai][Ai] += pcEqn;

    // alpha in pEqn [p, Ai] (already set)
    eqnBlock[0][Ai] += pcEqn;
}


// ************************************************************************* //
