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
    Foam::UEqnModels::Darcy

Description
    Segregated UEqn model for use in regions with no obvious dispersed phase.

    Reference:
    \verbatim

    \endverbatim

SourceFiles
    Darcy.C

\*---------------------------------------------------------------------------*/

#ifndef Darcy_H
#define Darcy_H

#include "UEqnModel.H"
#include "relativePhaseInterface.H"
#include "KrModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace UEqnModels
{

/*---------------------------------------------------------------------------*\
                         Class Darcy Declaration
\*---------------------------------------------------------------------------*/

class Darcy
:
    public UEqnModel
{
    // Private Data

        //- Interface
        const relativePhaseInterface interface_;

        //- Kr model
        autoPtr<KrModel> KrModel_;

        //- Darcy term without relative permeability
        surfaceScalarField Tauf;

        //- Expoent of saturation in relative permeability
        surfaceScalarField etaf_;


public:

    //- Runtime type information
    TypeName("Darcy");


    // Constructors

        //- Construct from a dictionary and an interface
        Darcy
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~Darcy();


    // Member Functions

        //- Access the interface
        virtual const phaseInterface& interface() const;

        //- Correct the model
        virtual void correct();

        //- Update the model after solution
        virtual void update();

        //- Insert the Taus fields
        virtual void insertTaus
        (
            HashPtrTable<surfaceScalarField>& taus
        );

        //- Add the matrices for the coupled algorithm
        virtual void addTerms
        (
            phaseSystem::eqnBlockTable& eqnBlock,
            PtrList<surfaceScalarField>& phi0s
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace UEqnModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
