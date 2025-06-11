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
    Foam::UEqnModel

Description
    Model for UEqn between phases

SourceFiles
    UEqnModel.C
    UEqnModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef UEqnModel_H
#define UEqnModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "phaseSystem.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class UEqnModel Declaration
\*---------------------------------------------------------------------------*/

class UEqnModel
:
    public regIOobject
{
public:

    //- Runtime type information
    TypeName("UEqnModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            UEqnModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseInterface& interface
            ),
            (dict, interface)
        );


    // Constructors

        // Construct from a dictionary and an interface
        UEqnModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~UEqnModel();


    // Selectors

        static autoPtr<UEqnModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Member Functions

        //- Access the interface
        virtual const phaseInterface& interface() const = 0;

        //- Correct the momentum
        virtual void correct() = 0;

        //- Update the model
        virtual void update() = 0;

        //- Insert the Taus fields
        virtual void insertTaus
        (
            HashPtrTable<surfaceScalarField>& taus
        ) = 0;

        //- Add the matrices for the coupled algorithm
        virtual void addTerms
        (
            phaseSystem::eqnBlockTable& eqnBlock,
            PtrList<surfaceScalarField>& phiHbyAs
        ) = 0;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
