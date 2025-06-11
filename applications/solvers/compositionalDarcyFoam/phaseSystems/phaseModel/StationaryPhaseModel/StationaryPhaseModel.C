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

#include "StationaryPhaseModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::StationaryPhaseModel<BasePhaseModel>::StationaryPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index)
{
    this->alphafRef() = fvc::interpolate(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::StationaryPhaseModel<BasePhaseModel>::~StationaryPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
bool Foam::StationaryPhaseModel<BasePhaseModel>::stationary() const
{
    return true;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::StationaryPhaseModel<BasePhaseModel>::U() const
{
    return volVectorField::null();
}


template<class BasePhaseModel>
Foam::volVectorField&
Foam::StationaryPhaseModel<BasePhaseModel>::URef()
{
    FatalErrorInFunction
        << "Cannot access the velocity of a stationary phase"
        << exit(FatalError);

    return const_cast<volVectorField&>(volVectorField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::phi() const
{
    return surfaceScalarField::null();
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StationaryPhaseModel<BasePhaseModel>::phiRef()
{
    FatalErrorInFunction
        << "Cannot access the flux of a stationary phase"
        << exit(FatalError);

    return const_cast<surfaceScalarField&>(surfaceScalarField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return surfaceScalarField::null();
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StationaryPhaseModel<BasePhaseModel>::alphaPhiRef()
{
    FatalErrorInFunction
        << "Cannot access the volumetric flux of a stationary phase"
        << exit(FatalError);

    return const_cast<surfaceScalarField&>(surfaceScalarField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::alphaRhoPhi() const
{
    return surfaceScalarField::null();
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::StationaryPhaseModel<BasePhaseModel>::alphaRhoPhiRef()
{
    FatalErrorInFunction
        << "Cannot access the mass flux of a stationary phase"
        << exit(FatalError);

    return const_cast<surfaceScalarField&>(surfaceScalarField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::StationaryPhaseModel<BasePhaseModel>::mu() const
{
    return volScalarField::null();
}

// ************************************************************************* //
