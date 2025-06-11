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
    Foam::MassTransferPhaseSystem

Description


SourceFiles
    MassTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef MassTransferPhaseSystem_H
#define MassTransferPhaseSystem_H

#include "phaseSystem.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class massTransferModel;

/*---------------------------------------------------------------------------*\
                 Class MassTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class MassTransferPhaseSystem
:
    public BasePhaseSystem
{
//private:

    // Private typedefs

        typedef HashTable
        <
            autoPtr<massTransferModel>,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > dmdtModelTable;

    // Private Data

        //- rAU field at faces (scalar)
        //HashPtrTable<surfaceScalarField> taus_;

        // Sub Models

            //- Mass transfer models
            dmdtModelTable dmdtModels_;


public:

    // Constructors

        //- Construct from fvMesh
        MassTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~MassTransferPhaseSystem();


    // Member Functions

        // Transfers

            //- Return the phiAlpha0P0s
            //virtual PtrList<surfaceScalarField>& phiAlpha0P0s();

            //- Add the matrices for the coupled algorithm
            virtual void addTerms
            (
                phaseSystem::eqnBlockTable& eqnBlock
            );


        // Evolution

            //- Correct the mass transfer models
            virtual void correct();

            //- Update the mass transfer models after solution
            virtual void update();


        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MassTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
