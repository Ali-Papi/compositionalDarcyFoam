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

#include "NonPorousPhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::NonPorousPhaseModel<BasePhaseModel>::NonPorousPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    S_
    (
        IOobject
        (
            IOobject::groupName("S", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::NonPorousPhaseModel<BasePhaseModel>::~NonPorousPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
bool Foam::NonPorousPhaseModel<BasePhaseModel>::porous() const
{
    return false;
}


template<class BasePhaseModel>
Foam::volScalarField&
Foam::NonPorousPhaseModel<BasePhaseModel>::SRef()
{
    return S_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::NonPorousPhaseModel<BasePhaseModel>::S() const
{
    return S_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volDiagTensorField>
Foam::NonPorousPhaseModel<BasePhaseModel>::invK() const
{
    return volDiagTensorField::null();
}


template<class BasePhaseModel>
Foam::volDiagTensorField&
Foam::NonPorousPhaseModel<BasePhaseModel>::invKRef()
{
    FatalErrorInFunction
        << "Cannot access the invK of a non-porous phase"
        << exit(FatalError);

    return const_cast<volDiagTensorField&>(volDiagTensorField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceDiagTensorField>
Foam::NonPorousPhaseModel<BasePhaseModel>::invKTf() const
{
    return surfaceDiagTensorField::null();
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::NonPorousPhaseModel<BasePhaseModel>::invKf() const
{
    return surfaceScalarField::null();
}

/*
template<class BasePhaseModel>
void Foam::NonPorousPhaseModel<BasePhaseModel>::rockThermo()
{}


template<class BasePhaseModel>
void Foam::NonPorousPhaseModel<BasePhaseModel>::correctInitialAlpha()
{}
*/
// ************************************************************************* //
