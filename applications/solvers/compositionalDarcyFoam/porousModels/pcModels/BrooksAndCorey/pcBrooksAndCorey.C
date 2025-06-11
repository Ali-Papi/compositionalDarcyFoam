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

#include "pcBrooksAndCorey.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace pcModels
{
    defineTypeNameAndDebug(pcBrooksAndCorey, 0);
    addToRunTimeSelectionTable(pcModel, pcBrooksAndCorey, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pcModels::pcBrooksAndCorey::pcBrooksAndCorey
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
    pc0(dict.lookupOrDefault("pc0", 0.)),
    zeta(dict.lookupOrDefault("zeta", 0.)),
    eta(dict.lookupOrDefault("eta", 0.))
{
    if (!dict.found("pc0") || !dict.found("eta") || !dict.found("zeta"))
    {
        FatalErrorInFunction
            << "BrooksAndCorey: pc0, zeta or eta not found."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pcModels::pcBrooksAndCorey::~pcBrooksAndCorey()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseInterface&
Foam::pcModels::pcBrooksAndCorey::interface() const
{
    return interface_;
}


void Foam::pcModels::pcBrooksAndCorey::correct()
{
    const phaseModel& phase = interface_.relative();
    const phaseSystem& fluid = interface_.fluid();
    const volScalarField& alpha = phase;
    const volScalarField& alphaVoid = fluid.alphaVoid();

    surfaceScalarField Sf(fvc::interpolate(alpha/alphaVoid));

    scalar iota = eta - zeta -1;

    forAll(KrdPcdSf, facei)
    {
        KrdPcdSf[facei] = -zeta*pc0*pow(max(Sf[facei], 0.001), iota); 
    }

    surfaceScalarField::GeometricBoundaryField& KrdPcdSfb = KrdPcdSf.boundaryField();
    surfaceScalarField::GeometricBoundaryField& Sfb = Sf.boundaryField();

    forAll(KrdPcdSfb, patchi)
    {
        if (!KrdPcdSfb[patchi].coupled()) continue;

        forAll(KrdPcdSfb[patchi], facei)
        {
            KrdPcdSfb[patchi][facei] = -zeta*pc0*pow(max(Sfb[patchi][facei], 0.001), iota);
        }
    }
}


void Foam::pcModels::pcBrooksAndCorey::update()
{}


void Foam::pcModels::pcBrooksAndCorey::addTerms
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
