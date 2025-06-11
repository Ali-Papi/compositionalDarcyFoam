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

#include "PorousPhaseModel.H"
#include "phaseSystem.H"
#include "reverseLinear.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::PorousPhaseModel<BasePhaseModel>::PorousPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    invK_
    (
        IOobject
        (
            IOobject::groupName("invK", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    invKTf_
    (
        IOobject
        (
            IOobject::groupName("invKTf", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        reverseLinear<diagTensor>(fluid.mesh()).interpolate(invK_)
    ),
    invKf_
    (
        IOobject
        (
            IOobject::groupName("invKf", this->name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (
            mag
            (
                (fluid.mesh().Sf()/mag(fluid.mesh().Sf()) & invKTf_) 
              & (fluid.mesh().Sf()/mag(fluid.mesh().Sf()))
            )
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::PorousPhaseModel<BasePhaseModel>::~PorousPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
bool Foam::PorousPhaseModel<BasePhaseModel>::porous() const
{
    return true;
}


template<class BasePhaseModel>
Foam::volScalarField&
Foam::PorousPhaseModel<BasePhaseModel>::SRef()
{
    FatalErrorInFunction
        << "Cannot access the S of a porous phase"
        << exit(FatalError);

    return const_cast<volScalarField&>(volScalarField::null());
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::PorousPhaseModel<BasePhaseModel>::S() const
{
    return volScalarField::null();
}


template<class BasePhaseModel>
Foam::tmp<Foam::volDiagTensorField>
Foam::PorousPhaseModel<BasePhaseModel>::invK() const
{
    return invK_;
}


template<class BasePhaseModel>
Foam::volDiagTensorField&
Foam::PorousPhaseModel<BasePhaseModel>::invKRef()
{
    return invK_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceDiagTensorField>
Foam::PorousPhaseModel<BasePhaseModel>::invKTf() const
{
    return invKTf_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::PorousPhaseModel<BasePhaseModel>::invKf() const
{
    return invKf_;
}

/*
template<class BasePhaseModel>
void Foam::PorousPhaseModel<BasePhaseModel>::rockThermo()
{
    if (!this->incompressible())
    {
        volScalarField& psi = const_cast<volScalarField&>(this->thermoRef().psi());
        volScalarField& rho = const_cast<volScalarField&>(this->thermoRef().rho());
        
        const volScalarField& p_ =  this->thermo().p();

        rho = rhoRef_/(1 - cf_*(p_ - pRef_));
        
        psi = cf_*rhoRef_/pow((cf_*(pRef_ -p_) + 1), 2);
    }
}

template<class BasePhaseModel>
void Foam::PorousPhaseModel<BasePhaseModel>::correctInitialAlpha()
{
    if (!this->incompressible())
    {
        const volScalarField& rho = this->thermoRef().rho();        
        volScalarField& alpha = *this;

        alpha = alpha0_*rhoRef_/rho;
    }
}
*/

// ************************************************************************* //
