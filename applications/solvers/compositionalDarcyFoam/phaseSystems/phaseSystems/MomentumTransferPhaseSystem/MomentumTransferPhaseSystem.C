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

#include "MomentumTransferPhaseSystem.H"

#include "UEqnModel.H"
#include "pcModel.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::
MomentumTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh)
{
    this->generateInterfacialModels(this->subDict("momentum"), UEqnModels_);
    this->generateInterfacialModels(this->subDict("capillarity"), pcModels_);

    if (pcModels_.size() > 0)
    {
        this->setCapillarity(true);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::
~MomentumTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::correct()
{
    BasePhaseSystem::correct();

    taus_.clear();

    surfaceScalarField alphaTausp
    (
        IOobject
        (
            "alphaTausp",
            this->mesh().time().timeName(),
            this->mesh()
        ),
        this->mesh(),
        dimensionedScalar("zero", dimensionSet(-1, 3, 1, 0, 0), Zero)
    );

    taus_.set
    (
        "alphaTausp",
        new surfaceScalarField(alphaTausp)
    );

    // Inserting phiHbyAs
    forAll(this->movingPhases(), movingPhasei)
    {
        const phaseModel& phase = this->movingPhases()[movingPhasei];

        if (phase.pure())
        {
            label Ai = phase.blockIndex()[0];

            phi0s_.set
            (
                Ai,
                new surfaceScalarField
                (
                    IOobject
                    (
                        IOobject::groupName("phi0", phase.name()),
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh(),
                    dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0), Zero)
                )
            );
        }
        else
        {
            forAll(phase.Y(), compi)
            {
                label Ai = phase.blockIndex()[compi];

                phi0s_.set
                (
                    Ai,
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName("phi0", phase.species()[compi]),
                            this->mesh().time().timeName(),
                            this->mesh()
                        ),
                        this->mesh(),
                        dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0), Zero)
                    )
                );
            }
        }
    }

    // Correct the UEqn models
    forAllIter
    (
        UEqnModelTable,
        UEqnModels_,
        UEqnModelIter
    )
    {
        UEqnModelIter()->correct();
        UEqnModelIter()->insertTaus(taus_);
    }

    // Correct the capillarity models
    forAllIter
    (
        pcModelTable,
        pcModels_,
        pcModelIter
    )
    {
        pcModelIter()->correct();
    }
}


template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::update()
{
    BasePhaseSystem::update();

    // Update the UEqn models
    forAllIter
    (
        UEqnModelTable,
        UEqnModels_,
        UEqnModelIter
    )
    {
        UEqnModelIter()->update();
    }

    // Update the capillarity models
    forAllIter
    (
        pcModelTable,
        pcModels_,
        pcModelIter
    )
    {
        pcModelIter()->update();
    }
}


template<class BasePhaseSystem>
void Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::addTerms
(
    phaseSystem::eqnBlockTable& eqnBlock
)
{
    BasePhaseSystem::addTerms(eqnBlock);

    // Add the UEqn models terms
    forAllIter
    (
        UEqnModelTable,
        UEqnModels_,
        UEqnModelIter
    )
    {
        UEqnModelIter()->addTerms(eqnBlock, phi0s_);
    }
    
    // Add the pc models terms
    forAllIter
    (
        pcModelTable,
        pcModels_,
        pcModelIter
    )
    {
        pcModelIter()->addTerms(eqnBlock, phi0s_);
    }
}


template<class BasePhaseSystem>
const Foam::surfaceScalarField&
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::taus
(
    const word& fieldName
) const
{
    return *taus_[fieldName];
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::surfaceScalarField>&
Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::phi0s()
{
    return phi0s_;
}


template<class BasePhaseSystem>
bool Foam::MomentumTransferPhaseSystem<BasePhaseSystem>::read()
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
