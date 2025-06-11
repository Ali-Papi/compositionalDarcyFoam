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

#include "blockSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockSystem, 0);
    defineRunTimeSelectionTable(blockSystem, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockSystem::blockSystem
(
    const label N,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_(dict),
    converged_(false),
    finalIter_(false),
    iter_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockSystem::~blockSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::blockSystem::dict()
{
    return dict_;
}


bool Foam::blockSystem::loop()
{
    const dictionary& coupledDict = mesh_.solutionDict().subDict("COUPLED");
    label nIterCOUPLED_ = coupledDict.lookupOrDefault("nOuterCorrectors", 1);

    iter_++;
    
    if (iter_ == nIterCOUPLED_ + 1)
    {
        Info<<"COUPLED: not converged within "
            << nIterCOUPLED_ << " iterations" << endl;

        iter_ = 0;
        return false;
    }
    
    bool completed = false;
    if (converged_)
    {
        if (finalIter_)
        {
            Info<< "COUPLED: converged in " << iter_ - 1
                << " iterations" << endl;

            iter_ = 0;
            converged_ = false;

            completed = true;
        }
        else
        {
            Info<< "COUPLED: iteration " << iter_ << endl;

            converged_ = true;
        }
    }
    else
    {
        if (iter_ <= nIterCOUPLED_)
        {
            Info<< "COUPLED: iteration " << iter_ << endl;

            completed = false;
        }
    }

    return !completed;
}


void Foam::blockSystem::checkConvergence
(
    const scalar& maxRes
)
{
    const dictionary& coupledDict = mesh_.solutionDict().subDict("COUPLED");
    scalar residualControl = coupledDict.subDict("residualControl").subDict("pAs").lookupOrDefault("tolerance", 0.0);

    if (maxRes > residualControl)
    {
        converged_ = false;
        finalIter_ = false;
    }
    else
    {
        if (converged_)
        {
            finalIter_ = true;
        }
        else
        {
            converged_ = true;
            finalIter_ = false;
        }
    }
}





// ************************************************************************* //
