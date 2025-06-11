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

#include "IsothermalPhaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::IsothermalPhaseModel<BasePhaseModel>::IsothermalPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, referencePhase, index)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::IsothermalPhaseModel<BasePhaseModel>::~IsothermalPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::IsothermalPhaseModel<BasePhaseModel>::correctThermo()
{
    BasePhaseModel::correctThermo();

    // Correct the thermo, but make sure that the temperature remains the same
/*    tmp<volScalarField> TCopy
    (
        new volScalarField
        (
            this->thermo().T().name() + ":Copy",
            this->thermo().T()
        )
    );
    this->thermo_->he() = this->thermo().he(this->thermo().p(), TCopy);*/
//    this->thermo_->correct();
//    this->thermo_->T() = TCopy;
}


template<class BasePhaseModel>
bool Foam::IsothermalPhaseModel<BasePhaseModel>::isothermal() const
{
    return true;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::IsothermalPhaseModel<BasePhaseModel>::heEqn()
{
    FatalErrorInFunction
        << "Cannot construct an energy equation for an isothermal phase"
        << exit(FatalError);

    return tmp<fvScalarMatrix>();
}


// ************************************************************************* //
