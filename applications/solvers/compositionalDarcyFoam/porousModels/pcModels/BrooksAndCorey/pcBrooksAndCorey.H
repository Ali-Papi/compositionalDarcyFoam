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
    Foam::pcModels::pcBrooksAndCorey

Description
    Standard Brooks and Corey capillary pressure model.

    This model must be used in conjunction with the BrooksAndCorey
    relative permeability model. Where the 'eta' parameter is the 
    saturation expoent in the relative permeability model and the 'zeta'
    parameter is the saturation exponent in this model.

Usage
    Example of model specification:
    \verbatim
    phaseInterface
    {
        type            BrooksAndCorey;
        pc0             0.001;
        zeta            0.5;
        eta             2;
    }
    \endverbatim

SourceFiles
    pcBrooksAndCorey.C

\*---------------------------------------------------------------------------*/

#ifndef pcBrooksAndCorey_H
#define pcBrooksAndCorey_H

#include "pcModel.H"
#include "relativePhaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace pcModels
{

/*---------------------------------------------------------------------------*\
                         Class pcBrooksAndCorey Declaration
\*---------------------------------------------------------------------------*/

class pcBrooksAndCorey
:
    public pcModel
{
    // Private Data

        //- Interface
        const relativePhaseInterface interface_;

        //- SurfaceField for Kr*dPc/dS
        surfaceScalarField KrdPcdSf;

        //- Capillarity pressure at S = 1 
        scalar pc0;

        //- Exponent of the Brooks and Corey capillary pressure model
        scalar zeta;

        //- Exponent of the Brooks and Corey relative permeability model
        scalar eta;


public:

    //- Runtime type information
    TypeName("BrooksAndCorey");


    // Constructors

        //- Construct from a dictionary and an interface
        pcBrooksAndCorey
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~pcBrooksAndCorey();


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
