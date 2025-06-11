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

#include "KrTabulated.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace KrModels
{
    defineTypeNameAndDebug(KrTabulated, 0);
    addToRunTimeSelectionTable(KrModel, KrTabulated, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::KrModels::KrTabulated::KrTabulated
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    KrModel(dict, interface),
    Krf_
    (
        IOobject
        (
            IOobject::groupName("Krf", interface.phase1().name()),
            interface.phase1().mesh().time().timeName(),
            interface.phase1().mesh()
        ),
        interface.phase1().mesh(),
        dimensionedScalar("zero", dimless, 0)
    ),
    etaf_
    (
        IOobject
        (
            IOobject::groupName("etaf", interface.phase1().name()),
            interface.phase1().mesh().time().timeName(),
            interface.phase1().mesh()
        ),
        interface.phase1().mesh(),
        dimensionedScalar("zero", dimless, 0)
    ),
    Kr_(dict)
{
    const phaseModel& phase = interface_.phase1();
    const volScalarField& alpha = phase;
    const volScalarField& alphaVoid = phase.fluid().alphaVoid();

    surfaceScalarField Seff(fvc::interpolate(alpha/alphaVoid));

    forAll(Krf_, facei)
    {
        Krf_[facei] = Kr_(Seff[facei]);
    }

    surfaceScalarField::GeometricBoundaryField& Krfb = Krf_.boundaryField();
    surfaceScalarField::GeometricBoundaryField& Seffb = Seff.boundaryField();

    forAll(Krfb, patchi)
    {
        forAll(Krfb[patchi], facei)
        {
            Krfb[patchi][facei] = Kr_(Seffb[patchi][facei]);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::KrModels::KrTabulated::~KrTabulated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::KrModels::KrTabulated::Krf
(
    const surfaceScalarField& Sf
)
{
    forAll(Krf_, facei)
    {
        Krf_[facei] += 0.5*Kr_(Sf[facei]);
    }

    surfaceScalarField::GeometricBoundaryField& Krfb = Krf_.boundaryField();
    const surfaceScalarField::GeometricBoundaryField& Sfb = Sf.boundaryField();

    forAll(Krfb, patchi)
    {
        forAll(Krfb[patchi], facei)
        {
            Krfb[patchi][facei] += 0.5*Kr_(Sfb[patchi][facei]);
        }
    }

    Krf_ /= 1.5;

    return Krf_;
}


Foam::tmp<Foam::surfaceScalarField> Foam::KrModels::KrTabulated::etaf
(
    const surfaceScalarField& Sf
)
{
    forAll(etaf_, facei)
    {
        if (Sf[facei] > 0.01)
        {
            etaf_[facei] = Kr_.rateOfChange(Sf[facei])/Sf[facei];
        }
        else
        {
            etaf_[facei] = 1.;
        }
        
        if (etaf_[facei] < 1)
        {
            etaf_[facei] = 1;
        }
        else if (etaf_[facei] > 3)
        {
            etaf_[facei] = 3;
        }
    }

    return etaf_;
}


// ************************************************************************* //
