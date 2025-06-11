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

#include "KrBrooksAndCorey.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace KrModels
{
    defineTypeNameAndDebug(KrBrooksAndCorey, 0);
    addToRunTimeSelectionTable(KrModel, KrBrooksAndCorey, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::KrModels::KrBrooksAndCorey::KrBrooksAndCorey
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    KrModel(dict, interface),
    Krmax
    (
        "Krmax",
        dimless,
        dict.lookupOrDefault("Krmax", 1.0)
    ),
    eta_
    (
        "eta",
        dimless,
        dict.lookup("eta")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::KrModels::KrBrooksAndCorey::~KrBrooksAndCorey()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::KrModels::KrBrooksAndCorey::Krf
(
    const surfaceScalarField& Sf
)
{
    return pow(Sf, eta_)*Krmax;
}


Foam::tmp<Foam::surfaceScalarField> Foam::KrModels::KrBrooksAndCorey::etaf
(
    const surfaceScalarField& Sf
)
{
    const phaseModel& phase1 = interface_.phase1();

    tmp<surfaceScalarField> tetaf
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("etaf", phase1.name()),
                phase1.fluid().mesh().time().timeName(),
                phase1.fluid().mesh()
            ),
            phase1.fluid().mesh(),
            eta_
        )
    );

    return tetaf;
}


// ************************************************************************* //
