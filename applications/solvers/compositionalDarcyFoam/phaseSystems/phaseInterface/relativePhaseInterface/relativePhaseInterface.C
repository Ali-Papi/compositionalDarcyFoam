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

#include "relativePhaseInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    bool relativePhaseInterfaceAddedHeadSeparator =
        phaseInterface::addHeadSeparator(relativePhaseInterface::separator()[0]);

    bool relativePhaseInterfaceAddedOldSeparatorToSeparator =
        phaseInterface::addOldSeparatorToSeparator
        (
            "in",
            relativePhaseInterface::separator()[0]
        );
}

namespace Foam
{
    defineTypeNameWithName
    (
        relativePhaseInterface,
        separatorsToTypeName(separator()).c_str()
    );
    addToRunTimeSelectionTable(phaseInterface, relativePhaseInterface, word);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativePhaseInterface::relativePhaseInterface
(
    const phaseModel& relative,
    const phaseModel& reference
)
:
    phaseInterface(relative, reference),
    relative_(relative)
{}


Foam::relativePhaseInterface::relativePhaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(fluid, name),
    relative_(identifyPhases(fluid, name, separator()).first())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativePhaseInterface::~relativePhaseInterface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::relativePhaseInterface::name() const
{
    return relative().name() + '_' + separator()[0] + '_' + reference().name();
}


const Foam::phaseModel& Foam::relativePhaseInterface::relative() const
{
    return relative_;
}


const Foam::phaseModel& Foam::relativePhaseInterface::reference() const
{
    return otherPhase(relative_);
}

// ************************************************************************* //
