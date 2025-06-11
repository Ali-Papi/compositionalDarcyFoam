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
    Foam::relativePhaseInterface

Description
    Class to represent a interface between phases where the two phases are
    considered to be relative; that is, separated by a geometrically complex
    interface for which dispersed representations are inappropriate.

SourceFiles
    relativePhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef relativePhaseInterface_H
#define relativePhaseInterface_H

#include "phaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class relativePhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class relativePhaseInterface
:
    virtual public phaseInterface
{
private:

    // Private Data

        //- Relative phase
        const phaseModel& relative_;


public:

    //- Runtime type information
    TypeName("relativePhaseInterface");


    // Constructors

        //- Construct from phases
        relativePhaseInterface
        (
            const phaseModel& relative,
            const phaseModel& reference
        );

        //- Construct from fluid and name
        relativePhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~relativePhaseInterface();


    // Static Member Functions

        //- Return the separator that delimits this interface's name
        static wordList separator()
        {
            return wordList(1, "relativeWith");
        }


    // Member Functions

        //- Name
        virtual word name() const;


        // Access

            //- Relative phase
            const phaseModel& relative() const;

            //- Reference phase
            const phaseModel& reference() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
