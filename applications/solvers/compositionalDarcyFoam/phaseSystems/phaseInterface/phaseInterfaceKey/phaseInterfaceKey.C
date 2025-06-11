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

#include "phaseModel.H"
#include "phaseInterface.H"
#include "phaseInterfaceKey.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseInterfaceKey::hash::hash()
{}


Foam::phaseInterfaceKey::phaseInterfaceKey()
{}


Foam::phaseInterfaceKey::phaseInterfaceKey
(
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    Pair<word>
    (
        phaseInterface::getPhase1(phase1, phase2).name(),
        phaseInterface::getPhase2(phase1, phase2).name()
    )
{}


Foam::phaseInterfaceKey::phaseInterfaceKey(const phaseInterface& interface)
:
    phaseInterfaceKey(interface.phase1(), interface.phase2())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseInterfaceKey::~phaseInterfaceKey()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

Foam::label Foam::phaseInterfaceKey::hash::operator()
(
    const phaseInterfaceKey& key
) const
{
    return word::hash()(key.first()) + word::hash()(key.second());
}


// * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * * //

bool Foam::operator==
(
    const phaseInterfaceKey& a,
    const phaseInterfaceKey& b
)
{
    return Pair<word>::compare(a, b) == 1;
}


bool Foam::operator!=
(
    const phaseInterfaceKey& a,
    const phaseInterfaceKey& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, phaseInterfaceKey& key)
{
    NotImplemented;

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const phaseInterfaceKey& key)
{
    NotImplemented;

    return os;
}


// ************************************************************************* //
