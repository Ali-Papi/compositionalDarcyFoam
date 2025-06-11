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
    Foam::NonPorousPhaseModel

Description
    Class which represents a non-porous phase, with no porousModels.

SourceFiles
    NonPorousPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef NonPorousPhaseModel_H
#define NonPorousPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NonPorousPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class NonPorousPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Saturation field
        volScalarField S_;
        
        
public:

    // Constructors

        NonPorousPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~NonPorousPhaseModel();


    // Member Functions

        // Momentum

            //- Return whether the phase is porous
            virtual bool porous() const;

            //- Return the saturation
            virtual tmp<volScalarField> S() const;

            //- Access the saturation
            virtual volScalarField& SRef();
            
            //- Return the absolut permeability
            virtual tmp<volDiagTensorField> invK() const;

            //- Access the absolut permeability
            virtual volDiagTensorField& invKRef();

            //- Return the absolut permeability vectorField at faces
            virtual tmp<surfaceDiagTensorField> invKTf() const;

            //- Return the absolut permeability at faces
            virtual tmp<surfaceScalarField> invKf() const;

            //- Correct the rock thermodynamics
            //virtual void rockThermo();

            //- Correct the volume fraction at initial time
            //virtual void correctInitialAlpha();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NonPorousPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
