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
    Foam::blockSystems::block3

Description

    Reference:
    \verbatim

    \endverbatim

SourceFiles
    block3.C

\*---------------------------------------------------------------------------*/

#ifndef block3_H
#define block3_H

#include "blockSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockSystems
{

/*---------------------------------------------------------------------------*\
                         Class block3 Declaration
\*---------------------------------------------------------------------------*/

class block3
:
    public blockSystem
{
    // Private Data

        //- Block field
        volVector3Field pAs_;
        
        //- Block matrix
        autoPtr<fvBlockMatrix<vector3>> pAsEqn_;


public:

    //- Runtime type information
    TypeName("block3");


    // Constructors

        //- Construct from a dictionary and an interface
        block3
        (
            const label N,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~block3();


    // Member Functions


        //- Create a new fvBlockMatrix
        virtual void newMatrix();


        // Insertion and retrieval public tools

            //- Retrieve part of internal field from this fvBlockMatrix
            virtual void retrieveSolution
            (
                const direction dir,
                scalarField& xSingle
            ) const;

            //- Insert matrix into this fvBlockMatrix
            virtual void insertEquation
            (
                const direction dir,
                fvScalarMatrix& matrix
            );

            //- Insert scalar equation coupling into this fvBlockMatrix
            virtual void insertEquationCoupling
            (
                const direction dirI,
                const direction dirJ,
                const fvScalarMatrix& matrix
            );

            //- Add field into block field
            virtual void blockAdd
            (
                const direction dir,
                const scalarField& xSingle
            );


        // Solver calls for fvBlockMatrix

            //- Solver controls read from fvSolution
            virtual void solve();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockSystems
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
