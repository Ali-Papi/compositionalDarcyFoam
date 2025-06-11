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

#include "KrMultiCompRockBrooksAndCorey.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace KrModels
{
    defineTypeNameAndDebug(KrMultiCompRockBrooksAndCorey, 0);
    addToRunTimeSelectionTable(KrModel, KrMultiCompRockBrooksAndCorey, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::KrModels::KrMultiCompRockBrooksAndCorey::KrMultiCompRockBrooksAndCorey
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    KrModel(dict, interface)
{
    const phaseModel& phase2 = interface.phase2();
 
    if (phase2.porous() && !phase2.pure())
    {
        forAll(phase2.Y(), compi)
        {
            const dictionary& compDict = dict.subDict(phase2.species()[compi]);

            eta_.append(compDict.lookupOrDefault("eta", 2));
            Krmax_.append(compDict.lookupOrDefault("Krmax", 1.0));
        }
    }
    else
    {
        FatalErrorInFunction
            << "Reference phase is not porous and multiComponent." << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::KrModels::KrMultiCompRockBrooksAndCorey::~KrMultiCompRockBrooksAndCorey()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::KrModels::KrMultiCompRockBrooksAndCorey::Krf
(
    const surfaceScalarField& Sf
)
{
    const phaseModel& phase1 = interface_.phase1();
    const phaseModel& phase2 = interface_.phase2();
    const volScalarField& alpha2 = phase2;

    volScalarField alpha2res(max(alpha2, phase2.residualAlpha()));

    tmp<surfaceScalarField> tKrf
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("Krf", phase1.name()),
                phase1.fluid().mesh().time().timeName(),
                phase1.fluid().mesh()
            ),
            phase1.fluid().mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    forAll(phase2.Y(), compi)
    {
        surfaceScalarField Yif(fvc::interpolate(phase2.Y()[compi]/alpha2res));

        tKrf() += Yif*pow(Sf, eta_[compi])*Krmax_[compi];    
    }

    return tKrf;
}


Foam::tmp<Foam::surfaceScalarField> Foam::KrModels::KrMultiCompRockBrooksAndCorey::etaf
(
    const surfaceScalarField& Sf
)
{
    const phaseModel& phase1 = interface_.phase1();
    const phaseModel& phase2 = interface_.phase2();
    const volScalarField& alpha2 = phase2;

    volScalarField alpha2res(max(alpha2, phase2.residualAlpha()));

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
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    forAll(phase2.Y(), compi)
    {
        surfaceScalarField Yif(fvc::interpolate(phase2.Y()[compi]/alpha2res));

        tetaf() += Yif*eta_[compi];    
    }

    return tetaf;
}


// ************************************************************************* //
