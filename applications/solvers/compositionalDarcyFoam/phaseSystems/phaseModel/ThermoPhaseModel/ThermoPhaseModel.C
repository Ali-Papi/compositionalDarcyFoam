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

#include "ThermoPhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel, class ThermoModel>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::ThermoPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index),
    thermo_(ThermoModel::New(fluid.mesh()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel, class ThermoModel>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::~ThermoPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class BasePhaseModel, class ThermoModel>
const Foam::basicRhoThermo&
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::thermo() const
{
    return thermo_();
}


template<class BasePhaseModel, class ThermoModel>
Foam::basicRhoThermo&
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::thermoRef()
{
    return thermo_();
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::volScalarField>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::rho() const
{
    return thermo_->rho();
}


template<class BasePhaseModel, class ThermoModel>
Foam::tmp<Foam::volScalarField>
Foam::ThermoPhaseModel<BasePhaseModel, ThermoModel>::mu() const
{
    return thermo_->mu();
}


// ************************************************************************* //
