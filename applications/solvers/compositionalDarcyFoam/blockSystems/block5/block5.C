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

\*---------------------------------------------------------------------------*/

#include "block5.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blockSystems
{
    defineTypeNameAndDebug(block5, 0);
    addToRunTimeSelectionTable(blockSystem, block5, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockSystems::block5::block5
(
    const label N,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    blockSystem(N, dict, mesh),
    pAs_
    (
        IOobject
        (
            "pAs",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector5("zero", dimless, vector5::zero)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockSystems::block5::~block5()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockSystems::block5::newMatrix()
{
    pAsEqn_.reset(new fvBlockMatrix<vector5>(pAs_));
}


void Foam::blockSystems::block5::retrieveSolution
(
    const direction dir,
    scalarField& xSingle
) const
{
    pAsEqn_->retrieveSolution(dir, xSingle);
}


void Foam::blockSystems::block5::insertEquation
(
    const direction dir,
    fvScalarMatrix& matrix
)
{
    pAsEqn_->insertEquation(dir, matrix);
}


void Foam::blockSystems::block5::insertEquationCoupling
(
    const direction dirI,
    const direction dirJ,
    const fvScalarMatrix& matrix
)
{
    pAsEqn_->insertEquationCoupling(dirI, dirJ, matrix);
}


void Foam::blockSystems::block5::blockAdd
(
    const direction dir,
    const scalarField& xSingle
)
{
    forAll (xSingle, cellI)
    {
        pAs_[cellI](dir) += xSingle[cellI];
    }
}


void Foam::blockSystems::block5::solve()
{
    const BlockSolverPerformance<vector5>& residual = pAsEqn_->solve(this->dict());
    scalar maxResidual = cmptMax(residual.initialResidual());
    
    checkConvergence(maxResidual);
}


// ************************************************************************* //
