/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "phaseModel.H"
//#include "ThermoPhaseModel.H"
#include "IsothermalPhaseModel.H"
//#include "AnisothermalPhaseModel.H"
#include "PurePhaseModel.H"
#include "MultiComponentPhaseModel.H"
#include "InertPhaseModel.H"
//#include "ReactingPhaseModel.H"
#include "MovingPhaseModel.H"
#include "StationaryPhaseModel.H"
#include "PorousPhaseModel.H"
#include "NonPorousPhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        NonPorousPhaseModel
        <
            IsothermalPhaseModel
            <
                PurePhaseModel
                <
                    InertPhaseModel
                    <
                        MovingPhaseModel
                        <
                            //ThermoPhaseModel<phaseModel, basicRhoThermo>
                            phaseModel
                        >
                    >
                >
            >
        >
        pureIsothermalPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureIsothermalPhaseModel,
        phaseSystem,
        pureIsothermalPhaseModel
    );

    typedef
        PorousPhaseModel
        <
            IsothermalPhaseModel
            <
                PurePhaseModel
                <
                    InertPhaseModel
                    <
                        StationaryPhaseModel
                        <
                            //ThermoPhaseModel<phaseModel, basicRhoThermo>
                            phaseModel
                        >
                    >
                >
            >
        >
        purePorousIsothermalPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        purePorousIsothermalPhaseModel,
        phaseSystem,
        purePorousIsothermalPhaseModel
    );

    typedef
        NonPorousPhaseModel
        <
            IsothermalPhaseModel
            <
                MultiComponentPhaseModel
                <
                    InertPhaseModel
                    <
                        MovingPhaseModel
                        <
                            //ThermoPhaseModel<phaseModel, rhoReactionThermo>
                            phaseModel
                        >
                    >
                >
            >
        >
        multiComponentIsothermalPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multiComponentIsothermalPhaseModel,
        phaseSystem,
        multiComponentIsothermalPhaseModel
    );

    typedef
        PorousPhaseModel
        <
            IsothermalPhaseModel
            <
                MultiComponentPhaseModel
                <
                    InertPhaseModel
                    <
                        StationaryPhaseModel
                        <
                            //ThermoPhaseModel<phaseModel, basicRhoThermo>
                            phaseModel
                        >
                    >
                >
            >
        >
        multiComponentPorousIsothermalPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multiComponentPorousIsothermalPhaseModel,
        phaseSystem,
        multiComponentPorousIsothermalPhaseModel
    );
}

// ************************************************************************* //
