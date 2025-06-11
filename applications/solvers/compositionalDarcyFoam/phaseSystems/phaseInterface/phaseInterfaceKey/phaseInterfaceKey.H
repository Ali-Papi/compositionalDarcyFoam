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

Class
    Foam::phaseInterfaceKey

Description
    Word-pair based class used for keying interface models in hash tables.

SourceFiles
    phaseInterfaceKey.C

\*---------------------------------------------------------------------------*/

#ifndef phaseInterfaceKey_H
#define phaseInterfaceKey_H

#include "Pair.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class phaseModel;
class phaseInterface;
class phaseInterfaceKey;

bool operator==(const phaseInterfaceKey&, const phaseInterfaceKey&);
bool operator!=(const phaseInterfaceKey&, const phaseInterfaceKey&);

Istream& operator>>(Istream&, phaseInterfaceKey&);
Ostream& operator<<(Ostream&, const phaseInterfaceKey&);

/*---------------------------------------------------------------------------*\
                       Class phaseInterfaceKey Declaration
\*---------------------------------------------------------------------------*/

class phaseInterfaceKey
:
    public Pair<word>
{
public:

    // Public Classes

        //- Hashing class
        class hash
        :
            public Hash<phaseInterfaceKey>
        {
        public:

            // Constructors

                // Construct null
                hash();


            // Member Operators

                // Generate a hash from a phase interface key
                label operator()(const phaseInterfaceKey& key) const;
        };


public:

    // Constructors

        //- Construct null
        phaseInterfaceKey();

        //- Construct from phases
        phaseInterfaceKey(const phaseModel& phase1, const phaseModel& phase2);

        //- Construct from phase interface
        phaseInterfaceKey(const phaseInterface& interface);


    // Destructor
    virtual ~phaseInterfaceKey();


    // Friend Operators

        //- Test if keys are equal
        friend bool operator==
        (
            const phaseInterfaceKey& a,
            const phaseInterfaceKey& b
        );

        //- Test if keys are unequal
        friend bool operator!=
        (
            const phaseInterfaceKey& a,
            const phaseInterfaceKey& b
        );

        //- Read from stream
        friend Istream& operator>>(Istream& is, phaseInterfaceKey& key);

        //- Write to stream
        friend Ostream& operator<<(Ostream& os, const phaseInterfaceKey& key);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
