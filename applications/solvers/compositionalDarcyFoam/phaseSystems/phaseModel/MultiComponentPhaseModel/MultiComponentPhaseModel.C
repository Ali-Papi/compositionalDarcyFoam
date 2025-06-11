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

#include "MultiComponentPhaseModel.H"

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MultiComponentPhaseModel<BasePhaseModel>::MultiComponentPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -3, 0, 0, 0), 0)
    ),
    psi_
    (
        IOobject
        (
            IOobject::groupName("psi", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimensionSet(0, -2, 2, 0, 0), 0)
    ),
    species_(fluid.subDict(phaseName).lookup("species"))
{
    Y_.setSize(species_.size());
    psis_.setSize(species_.size());
    rho0s_.setSize(species_.size());
    rhos_.setSize(species_.size());

    forAll(Y_, i)
    {
        Y_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(species_[i], phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluid.mesh()
            )
        );

        psis_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("psi", species_[i]+"_In_"+phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar
                (
                    "zero",
                    dimensionSet(0, -2, 2, 0, 0),
                    fluid.subDict(phaseName).lookup("psi."+species_[i])
                )
            )
        );

        rho0s_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("rho0", species_[i]+"_In_"+phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid.mesh(),
                dimensionedScalar
                (
                    "zero",
                    dimensionSet(1, -3, 0, 0, 0),
                    fluid.subDict(phaseName).lookup("rho0."+species_[i])
                )
            )
        );

        rhos_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("rho", species_[i]+"_In_"+phaseName),
                    fluid.mesh().time().timeName(),
                    fluid.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                (psis_[i]*fluid.p() + rho0s_[i])
            )
        );
    }

    updateThermo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MultiComponentPhaseModel<BasePhaseModel>::~MultiComponentPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MultiComponentPhaseModel<BasePhaseModel>::updateThermo()
{
    const volScalarField& alpha = *this;

    forAll(rho_, celli)
    {
        if (alpha[celli] > 1e-10)
        {
            rho_[celli] = 0;
            psi_[celli] = 0;

            forAll(Y_, compi)
            {
                rho_[celli] += Y_[compi][celli]*rhos_[compi][celli];
                psi_[celli] += Y_[compi][celli]*psis_[compi][celli];
            }

            rho_[celli] /= alpha[celli];
            psi_[celli] /= alpha[celli];
        }
        else
        {
            rho_[celli] = rhos_[0][celli];
            psi_[celli] = psis_[0][celli];
        }
    }

    volScalarField::GeometricBoundaryField& rhob = rho_.boundaryField();
    volScalarField::GeometricBoundaryField& psib = psi_.boundaryField();
    const volScalarField::GeometricBoundaryField& alphab = alpha.boundaryField();

    forAll(rhob, patchi)
    {
        forAll(rhob[patchi], facei)
        {
            if (alphab[patchi][facei] > 1e-10)
            {
                rhob[patchi][facei] = 0;
                psib[patchi][facei] = 0;

                forAll(Y_, compi)
                {
                    rhob[patchi][facei] += 
                        Y_[compi].boundaryField()[patchi][facei]
                      * rhos_[compi].boundaryField()[patchi][facei];

                    psib[patchi][facei] += 
                        Y_[compi].boundaryField()[patchi][facei]
                      * psis_[compi].boundaryField()[patchi][facei];
                }

                rhob[patchi][facei] /= alphab[patchi][facei];
                psib[patchi][facei] /= alphab[patchi][facei];
            }
            else
            {
                rhob[patchi][facei] = rhos_[0].boundaryField()[patchi][facei];
                psib[patchi][facei] = psis_[0].boundaryField()[patchi][facei];
            } 
        }
    }
}


template<class BasePhaseModel>
void Foam::MultiComponentPhaseModel<BasePhaseModel>::correct()
{
    // update equation of state
    forAll(Y_, compi)
    {
        rhos_[compi] = psis_[compi]*this->fluid().p() + rho0s_[compi];
    }

    updateThermo();

    BasePhaseModel::correct();
}


template<class BasePhaseModel>
void Foam::MultiComponentPhaseModel<BasePhaseModel>::correctSpecies()
{
    BasePhaseModel::correctSpecies();
}


template<class BasePhaseModel>
bool Foam::MultiComponentPhaseModel<BasePhaseModel>::pure() const
{
    return false;
}


template<class BasePhaseModel>
const Foam::wordList&
Foam::MultiComponentPhaseModel<BasePhaseModel>::species() const
{
    return species_;
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::Y() const
{
    return Y_;
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::Yrho() const
{
    return rhos_;
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MultiComponentPhaseModel<BasePhaseModel>::Y(const word& name) const
{
    label Yi = 0;
    
    forAll(species_, speciei)
    {
        if (species_[speciei] == name)
        {
            Yi = speciei;
            break;
        }
    }

    return Y_[Yi];
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::YRef()
{
    return Y_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MultiComponentPhaseModel<BasePhaseModel>::rho() const
{
    return rho_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volScalarField>
Foam::MultiComponentPhaseModel<BasePhaseModel>::psi() const
{
    return psi_;
}


// ************************************************************************* //
