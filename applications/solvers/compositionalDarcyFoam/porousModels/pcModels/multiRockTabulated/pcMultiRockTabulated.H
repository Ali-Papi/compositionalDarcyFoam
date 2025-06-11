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
    Foam::pcModels::pcMultiRockTabulated

Description
    This is a mixed model of capillary pressure, i.e., for a multicomponent
    rock it is possible to specify a model for each component through a table.

    For each model, 2 tables are required, KrdPcdS and KrPcS tables (Kr*dpc/dS
    and Kr*pc/S, respectively). Relative permeability must be included in the
    table values.

Usage
    Example of model specification:
    \verbatim
    phaseInterface
    {
        type            multiRockTabulated;
        porousPhase     rock;
        KrdPcdS.specieI
        {
            fileName        "constant/KrdPcdSTable-I";
            outOfBounds     clamp;
        }
        KrPcS.specieI
        {
            fileName        "constant/KrPcSTable-I";
            outOfBounds     clamp;
        }
        ...
    }
    \endverbatim

SourceFiles
    pcMultiRockTabulated.C

\*---------------------------------------------------------------------------*/

#ifndef pcMultiRockTabulated_H
#define pcMultiRockTabulated_H

#include "pcModel.H"
#include "relativePhaseInterface.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pcModels
{

/*---------------------------------------------------------------------------*\
                         Class pcMultiRockTabulated Declaration
\*---------------------------------------------------------------------------*/

class pcMultiRockTabulated
:
    public pcModel
{
    // Private Data

        //- Interface
        const relativePhaseInterface interface_;

        //- Sum of Kr*dPc/dS of all rock components
        surfaceScalarField KrdPcdSf;

        //- Sum of Kr*grad(Yi)*Pc/S of all rock components
        surfaceScalarField KrgradYPcf;

        //- The Kr*dPc/dS table being used, including the bounding treatment
        PtrList<interpolationTable<scalar>> KrdPcdSTables;

        //- The Kr*Pc/S table being used, including the bounding treatment
        PtrList<interpolationTable<scalar>> KrPcSTables;

        //- Scale for adjustments 
        scalar scale_;

        //- Porous phase reference
        const phaseModel& porous_;


public:

    //- Runtime type information
    TypeName("multiRockTabulated");


    // Constructors

        //- Construct from a dictionary and an interface
        pcMultiRockTabulated
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~pcMultiRockTabulated();


    // Member Functions

        //- Access the interface
        virtual const phaseInterface& interface() const;

        //- Correct the model
        virtual void correct();

        //- Update the model after solution
        virtual void update();

        //- Add the matrices for the coupled algorithm
        virtual void addTerms
        (
            phaseSystem::eqnBlockTable& eqnBlock,
            PtrList<surfaceScalarField>& phi0s
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pcModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
