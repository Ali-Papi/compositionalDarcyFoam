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

#include "MassTransferPhaseSystem.H"

#include "massTransferModel.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::
MassTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generateInterfacialModels(this->subDict("massTransfer"), dmdtModels_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MassTransferPhaseSystem<BasePhaseSystem>::
~MassTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    // Correct mass transfer models
    forAllIter
    (
        dmdtModelTable,
        dmdtModels_,
        dmdtModelIter
    )
    {
        dmdtModelIter()->correct();
    }
}


template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::update()
{
    BasePhaseSystem::update();

    // Update mass transfer models
    forAllIter
    (
        dmdtModelTable,
        dmdtModels_,
        dmdtModelIter
    )
    {
        dmdtModelIter()->update();
    }
}


template<class BasePhaseSystem>
void Foam::MassTransferPhaseSystem<BasePhaseSystem>::addTerms
(
    phaseSystem::eqnBlockTable& eqnBlock
)
{
    BasePhaseSystem::addTerms(eqnBlock);

    // Adding the mass transfer models terms
    forAllIter
    (
        dmdtModelTable,
        dmdtModels_,
        dmdtModelIter
    )
    {
        dmdtModelIter()->addTerms(eqnBlock);
    }
}


template<class BasePhaseSystem>
bool Foam::MassTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Read models ...

        return readOK;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
