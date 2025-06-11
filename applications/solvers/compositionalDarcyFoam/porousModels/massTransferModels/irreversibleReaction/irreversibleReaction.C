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

#include "irreversibleReaction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferModels
{
    defineTypeNameAndDebug(irreversibleReaction, 0);
    addToRunTimeSelectionTable(massTransferModel, irreversibleReaction, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferModels::irreversibleReaction::irreversibleReaction
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    massTransferModel(dict, interface),
    interface_(interface.modelCast<massTransferModel, relativePhaseInterface>()),
    K_
    (
        "K",
        dimensionSet(1, -3, -1, 0, 0),
        dict.lookup("K")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massTransferModels::irreversibleReaction::~irreversibleReaction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseInterface&
Foam::massTransferModels::irreversibleReaction::interface() const
{
    return interface_;
}


void Foam::massTransferModels::irreversibleReaction::correct()
{}


void Foam::massTransferModels::irreversibleReaction::update()
{}


void Foam::massTransferModels::irreversibleReaction::addTerms
(
    phaseSystem::eqnBlockTable& eqnBlock
)
{
    const phaseModel& phase = interface_.relative();
    const volScalarField& alpha = phase;

    // Reaction: A -> B
    label Aa = phase.blockIndex()[0];
    label Ab = phase.blockIndex()[1];

    const volScalarField& Ya = phase.Y()[0];
    const volScalarField& rhoa = phase.Yrho()[0];
    const volScalarField& rhob = phase.Yrho()[1];

    eqnBlock[Aa][Aa] += fvm::Sp(K_/(alpha*rhoa),Ya);
    eqnBlock[0][Aa] += fvm::Sp(K_/(alpha*rhoa) -K_/(alpha*rhob), Ya);
    eqnBlock[Ab].set(Aa, fvm::Sp(-K_/(alpha*rhob),Ya));
}


// ************************************************************************* //
