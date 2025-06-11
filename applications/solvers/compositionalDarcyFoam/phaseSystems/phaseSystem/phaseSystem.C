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

#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseSystem, 0);
    defineRunTimeSelectionTable(phaseSystem, dictionary);
}

const Foam::word Foam::phaseSystem::propertiesName("phaseProperties");


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::phaseSystem::phiCalc() const
{
    tmp<surfaceScalarField> tmpPhi
    (
        new surfaceScalarField
        (
            IOobject
            (
                "tphi",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0), Zero)
        )
    );

    forAll(phaseModels_, movingPhasei)
    {
        const phaseModel& phase = phaseModels_[movingPhasei]; 

        if (phase.stationary()) continue;
 
        tmpPhi() += phase.alphaf()*phase.phi();
    }

    return tmpPhi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSystem::phaseSystem
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            propertiesName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    alphaVoid_
    (
        IOobject
        (
            "alphaVoid",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    ),
    porosity_
    (
        IOobject
        (
            "porosity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    ),
    snGradPSf_(nullptr),
    referencePhaseName_(lookupOrDefault("referencePhase", word::null)),
    phaseModels_
    (
        lookup("phases"),
        phaseModel::iNew(*this, referencePhaseName_)
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (phiCalc())
    ),
    phip_
    (
        IOobject
        (
            "phip",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    phig_
    (
        IOobject
        (
            "phig",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    mvConvP_(nullptr),
    mvConvG_(nullptr),
    hasCapillarity_(false)
{
    // Groupings
    label movingPhasei = 0;
    label stationaryPhasei = 0;
    label porousPhasei = 0;
    label anisothermalPhasei = 0;
    label multiComponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        movingPhasei += !phase.stationary();
        stationaryPhasei += phase.stationary();
        porousPhasei += phase.porous();
        anisothermalPhasei += !phase.isothermal();
        multiComponentPhasei += !phase.pure();
    }
    movingPhaseModels_.resize(movingPhasei);
    stationaryPhaseModels_.resize(stationaryPhasei);
    porousPhaseModels_.resize(porousPhasei);
    anisothermalPhaseModels_.resize(anisothermalPhasei);
    multiComponentPhaseModels_.resize(multiComponentPhasei);

    movingPhasei = 0;
    stationaryPhasei = 0;
    porousPhasei = 0;
    anisothermalPhasei = 0;
    multiComponentPhasei = 0;
    forAll(phaseModels_, phasei)
    {
        phaseModel& phase = phaseModels_[phasei];
        if (!phase.stationary())
        {
            movingPhaseModels_.set(movingPhasei++, &phase);
        }
        if (phase.stationary())
        {
            stationaryPhaseModels_.set(stationaryPhasei++, &phase);
        }
        if (phase.porous())
        {
            porousPhaseModels_.set(porousPhasei++, &phase);
        }
        if (!phase.isothermal())
        {
            anisothermalPhaseModels_.set(anisothermalPhasei++, &phase);
        }
        if (!phase.pure())
        {
            multiComponentPhaseModels_.set(multiComponentPhasei++, &phase);
        }
    }

    // Set the optional reference phase fraction from the other phases
    if (referencePhaseName_ != word::null)
    {
        phaseModel* referencePhasePtr = &phases()[0];
        volScalarField& referenceAlpha = *referencePhasePtr;

        referenceAlpha = 1;

        forAll(phaseModels_, phasei)
        {
            if (&phaseModels_[phasei] != referencePhasePtr)
            {
                referenceAlpha -= phaseModels_[phasei];
            }
        }
    }

    forAll(phases(), phasei)
    {
        const volScalarField& alphai = phases()[phasei];
        mesh_.schemesDict().setFluxRequired(alphai.name());
    }
    
    forAll(stationaryPhases(), phasei)
    {
        phaseModel& phase = stationaryPhases()[phasei];
        alphaVoid_ -= phase;
    }
    
    forAll(porousPhases(), phasei)
    {
        phaseModel& phase = porousPhases()[phasei];
        porosity_ -= phase;
    }

    forAll(phases(), phasei)
    {
        phaseModel& phase = phases()[phasei];
        
        if (phase.porous()) continue;
        
        const volScalarField& alphai = phase;
        phase.SRef() = alphai/porosity_;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseSystem::~phaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::phaseSystem::foundPhase
(
    const word& name
) const
{
    forAll(phaseModels_, i)
    {
        if (phaseModels_[i].name() == name)
        {
            return true;
        }
    }

    return false;
}


const Foam::phaseModel& Foam::phaseSystem::lookupPhase
(
    const word& name
) const
{
    forAll(phaseModels_, i)
    {
        if (phaseModels_[i].name() == name)
        {
            return phaseModels_[i];
        }
    }

    FatalErrorInFunction
        << "this phaseSystem does not contain phase " << name
        << exit(FatalError);

    return phaseModels_[0];
}


Foam::tmp<Foam::volScalarField> Foam::phaseSystem::rho() const
{
    tmp<volScalarField> trho(movingPhaseModels_[0]*movingPhaseModels_[0].rho());

    for
    (
        label movingPhasei=1;
        movingPhasei<movingPhaseModels_.size();
        movingPhasei++
    )
    {
        trho() +=
            movingPhaseModels_[movingPhasei]
           *movingPhaseModels_[movingPhasei].rho();
    }

    if (stationaryPhaseModels_.empty())
    {
        return trho;
    }
    else
    {
        return trho/alphaVoid_;
    }
}


Foam::tmp<Foam::volVectorField> Foam::phaseSystem::U() const
{
    tmp<volVectorField> tU(movingPhaseModels_[0]*movingPhaseModels_[0].U());

    for
    (
        label movingPhasei=1;
        movingPhasei<movingPhaseModels_.size();
        movingPhasei++
    )
    {
        tU() +=
            movingPhaseModels_[movingPhasei]
           *movingPhaseModels_[movingPhasei].U();
    }

    if (stationaryPhaseModels_.empty())
    {
        return tU;
    }
    else
    {
        return tU/alphaVoid_;
    }
}


void Foam::phaseSystem::correct()
{
    forAll(phaseModels_, phasei)
    {
        phaseModels_[phasei].correct();
    }
}


void Foam::phaseSystem::update()
{}


void Foam::phaseSystem::addTerms
(
    phaseSystem::eqnBlockTable& eqnBlock
)
{}


void Foam::phaseSystem::correctBoundaryFlux()
{
    forAll(movingPhases(), movingPhasei)
    {
        phaseModel& phase = movingPhases()[movingPhasei];

        const fvBoundaryMesh& boundary = mesh_.boundary();

        const volVectorField::GeometricBoundaryField& UBf = phase.U()().boundaryField();

        surfaceScalarField::GeometricBoundaryField& phiBf = phase.phiRef().boundaryField();

        forAll(boundary, patchi)
        {
            if (isA<fixedValueFvsPatchScalarField>(phiBf[patchi]))
            {
                phiBf[patchi] == (boundary[patchi].Sf() & UBf[patchi]);
            }
        }
    }
}

bool Foam::phaseSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        forAll(phaseModels_, phasei)
        {
            readOK &= phaseModels_[phasei].read();
        }

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}
// ************************************************************************* //
