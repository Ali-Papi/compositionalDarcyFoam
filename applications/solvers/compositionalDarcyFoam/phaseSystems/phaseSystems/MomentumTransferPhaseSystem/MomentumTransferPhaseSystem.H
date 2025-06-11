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
    Foam::MomentumTransferPhaseSystem

Description


SourceFiles
    MomentumTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef MomentumTransferPhaseSystem_H
#define MomentumTransferPhaseSystem_H

#include "phaseSystem.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class UEqnModel;
class pcModel;

/*---------------------------------------------------------------------------*\
                 Class MomentumTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class MomentumTransferPhaseSystem
:
    public BasePhaseSystem
{
private:

    // Private typedefs

        typedef HashTable
        <
            autoPtr<UEqnModel>,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > UEqnModelTable;

        typedef HashTable
        <
            autoPtr<pcModel>,
            phaseInterfaceKey,
            phaseInterfaceKey::hash
        > pcModelTable;

    // Private Data

        //- rAU field at faces (scalar)
        HashPtrTable<surfaceScalarField> taus_;

        //- phi0s fields (explicit fluxes)
        PtrList<surfaceScalarField> phi0s_;

        // Sub Models

            //- UEqn models
            UEqnModelTable UEqnModels_;

            //- pc models
            pcModelTable pcModels_;


public:

    // Constructors

        //- Construct from fvMesh
        MomentumTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~MomentumTransferPhaseSystem();


    // Member Functions

        // Transfers

            //- Return the Taus (rAUfs)
            virtual const surfaceScalarField& taus
            (
                const word& fieldName
            ) const;

            //- Return the phi0s
            virtual PtrList<surfaceScalarField>& phi0s();

            //- Add the matrices for the coupled algorithm
            virtual void addTerms
            (
                phaseSystem::eqnBlockTable& eqnBlock
            );


        // Evolution

            //- Correct the momemtum models
            virtual void correct();

            //- Update the momemtum models after solution
            virtual void update();


        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MomentumTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
