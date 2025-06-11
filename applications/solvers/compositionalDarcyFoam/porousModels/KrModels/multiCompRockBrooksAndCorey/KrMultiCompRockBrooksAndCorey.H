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
    Foam::KrModels::KrMultiCompRockBrooksAndCorey

Description
    Standard KrMultiCompRockBrooks and Corey relative permeability model.
    
    Reference:
    \verbatim

    \endverbatim

SourceFiles
    KrMultiCompRockBrooksAndCorey.C

\*---------------------------------------------------------------------------*/

#ifndef KrMultiCompRockBrooksAndCorey_H
#define KrMultiCompRockBrooksAndCorey_H

#include "KrModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace KrModels
{

/*---------------------------------------------------------------------------*\
                           Class KrMultiCompRockBrooksAndCorey Declaration
\*---------------------------------------------------------------------------*/

class KrMultiCompRockBrooksAndCorey
:
    public KrModel
{
protected:

    // Protected data

        //- Maximum relative permeability 
        scalarList Krmax_;

        //- Exponent of the BC model
        scalarList eta_;


public:

    //- Runtime type information
    TypeName("multiCompRockBrooksAndCorey");


    // Constructors

        //- Construct from a dictionary and a phase list
        KrMultiCompRockBrooksAndCorey
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~KrMultiCompRockBrooksAndCorey();


    // Member Functions

        //- Return the relative permeability
        virtual tmp<surfaceScalarField> Krf
        (
            const surfaceScalarField& Sf
        );

        //- Return the expoent of saturation in relative permeability
        virtual tmp<surfaceScalarField> etaf
        (
            const surfaceScalarField& Sf
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace KrModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
