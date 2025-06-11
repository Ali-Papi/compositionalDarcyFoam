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

#include "phaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseModel, 0);
    defineRunTimeSelectionTable(phaseModel, phaseSystem);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const bool referencePhase,
    const label index
)
:
    volScalarField
    (
        referencePhase
      ? volScalarField
        (
            IOobject
            (
                IOobject::groupName("alpha", phaseName),
                fluid.mesh().time().timeName(),
                fluid.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluid.mesh(),
            dimensionedScalar("zero", dimless, 0)
        )
      : volScalarField
        (
            IOobject
            (
                IOobject::groupName("alpha", phaseName),
                fluid.mesh().time().timeName(),
                fluid.mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluid.mesh()
        )
    ),
    fluid_(fluid),
    name_(phaseName),
    index_(index),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        fluid.subDict(phaseName).lookup("residualAlpha")
    ),
    alphaMax_(fluid.subDict(phaseName).lookupOrDefault("alphaMax", 1.0)),
    //divU_(nullptr),
    alphaf_
    (
        IOobject
        (
            IOobject::groupName("alphaf", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimless, 0) 
    )
{}


Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return autoPtr<phaseModel>(nullptr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::phaseModel::name() const
{
    return name_;
}


const Foam::word& Foam::phaseModel::keyword() const
{
    return name_;
}


Foam::label Foam::phaseModel::index() const
{
    return index_;
}


const Foam::phaseSystem& Foam::phaseModel::fluid() const
{
    return fluid_;
}


const Foam::dimensionedScalar& Foam::phaseModel::residualAlpha() const
{
    return residualAlpha_;
}


Foam::scalar Foam::phaseModel::alphaMax() const
{
    return alphaMax_;
}


Foam::surfaceScalarField& Foam::phaseModel::alphafRef()
{
    return alphaf_;
}


Foam::tmp<Foam::surfaceScalarField> Foam::phaseModel::alphaf() const
{
    return alphaf_;
}


void Foam::phaseModel::correct()
{}


void Foam::phaseModel::correctThermo()
{}


void Foam::phaseModel::correctReactions()
{}


void Foam::phaseModel::correctSpecies()
{}


bool Foam::phaseModel::read()
{
    return true;
}


const Foam::labelList& Foam::phaseModel::blockIndex() const
{
    return blockIndex_;
}


Foam::labelList& Foam::phaseModel::blockIndex()
{
    return blockIndex_;
}


void Foam::phaseModel::correctInflowOutflow(surfaceScalarField& alphaPhi) const
{
/*    surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();
    const volScalarField::Boundary& alphaBf = boundaryField();
    const surfaceScalarField::Boundary& phiBf = phi()().boundaryField();

    forAll(alphaPhiBf, patchi)
    {
        fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];

        if (!alphaPhip.coupled())
        {
            alphaPhip = phiBf[patchi]*alphaBf[patchi];
        }
    }*/
}

/*
Foam::tmp<Foam::volScalarField>
Foam::phaseModel::divU() const
{
    return divU_.valid() ? tmp<volScalarField>(divU_()) : tmp<volScalarField>();
}


void Foam::phaseModel::divU(tmp<volScalarField> divU)
{
    if (!divU_.valid())
    {
        divU_ = divU;
        divU_->rename(IOobject::groupName("divU", this->name()));
        divU_->checkIn();
    }
    else
    {
        divU_() = divU;
    }
}*/


// ************************************************************************* //
