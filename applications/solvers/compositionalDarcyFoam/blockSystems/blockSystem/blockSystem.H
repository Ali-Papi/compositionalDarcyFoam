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
    Foam::blockSystem

Description


SourceFiles
    blockSystem.C
    blockSystemNew.C

\*---------------------------------------------------------------------------*/

#ifndef blockSystem_H
#define blockSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class blockSystem Declaration
\*---------------------------------------------------------------------------*/

class blockSystem
{
protected:

    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the dict
        const dictionary& dict_;

        //- Convergence flag
        bool converged_;

        //- Final iteration flag
        bool finalIter_;

        //- Coupled iteration counter
        label iter_;


public:

    //- Runtime type information
    TypeName("blockSystem");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            blockSystem,
            dictionary,
            (
                const label N,
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (N, dict, mesh)
        );


    // Constructors

        // Construct from a dictionary and an interface
        blockSystem
        (
            const label N,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~blockSystem();


    // Selectors

        static autoPtr<blockSystem> New
        (
            const label N,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Return dictionary
        const dictionary& dict();

        //- Check convergence for coupled matrix
        void checkConvergence(const scalar& maxRes);

        //- COUPLED loop
        virtual bool loop();

        //- Create a new fvBlockMatrix
        virtual void newMatrix() = 0;


        // Insertion and retrieval public tools

            //- Retrieve part of internal field from this fvBlockMatrix
            virtual void retrieveSolution
            (
                const direction dir,
                scalarField& xSingle
            ) const = 0;

            //- Insert matrix into this fvBlockMatrix
            virtual void insertEquation
            (
                const direction dir,
                fvScalarMatrix& matrix
            ) = 0;

            //- Insert scalar equation coupling into this fvBlockMatrix
            virtual void insertEquationCoupling
            (
                const direction dirI,
                const direction dirJ,
                const fvScalarMatrix& matrix
            ) = 0;

            //- Add field into block field
            virtual void blockAdd
            (
                const direction dir,
                const scalarField& xSingle
            ) = 0;


        // Solver calls for fvBlockMatrix

            //- Solver controls read from fvSolution
            virtual void solve() = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
