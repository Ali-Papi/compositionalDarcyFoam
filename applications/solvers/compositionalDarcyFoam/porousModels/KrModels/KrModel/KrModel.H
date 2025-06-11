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
    Foam::KrModel

Description

SourceFiles
    KrModel.C
    newKrModel.C

\*---------------------------------------------------------------------------*/

#ifndef KrModel_H
#define KrModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseInterface.H"
#include "surfaceInterpolate.H"
#include "convectionScheme.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class KrModel Declaration
\*---------------------------------------------------------------------------*/

class KrModel
{
protected:

    // Protected data

        //- phase model reference        
        const phaseInterface& interface_;

public:

    //- Runtime type information
    TypeName("KrModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            KrModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface
            ),
            (dict, interface)
        );


    // Constructors

        //- Construct with dictionary name
        KrModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~KrModel();


    // Selectors

        //- Returns a reference to the selected Kr model
        static autoPtr<KrModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Return the relative permeability
        virtual tmp<surfaceScalarField> Krf
        (
            const surfaceScalarField& Sf
        ) = 0;

        //- Return the expoent of saturation in relative permeability
        virtual tmp<surfaceScalarField> etaf
        (
            const surfaceScalarField& Sf
        ) = 0;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
