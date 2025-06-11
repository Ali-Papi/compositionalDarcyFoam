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
    Foam::MultiComponentPhaseModel

Description
    Class which represents a phase with multiple species. Returns the species'
    mass fractions, and their governing equations.

SourceFiles
    MultiComponentPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef MultiComponentPhaseModel_H
#define MultiComponentPhaseModel_H

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class MultiComponentPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class MultiComponentPhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Density field
        volScalarField rho_;

        //- Compressibility field
        volScalarField psi_;

        //- List of species name
        wordList species_;

        //- Pointer list of species
        PtrList<volScalarField> Y_;

        //- Pointer list of densities
        PtrList<volScalarField> rhos_;

        //- Equation of state: rho = psi*p + rho0
        //- Linear coefficient of rho
        PtrList<volScalarField> psis_;

        //- Angular coefficient of rho
        PtrList<volScalarField> rho0s_;


public:

    // Constructors

        MultiComponentPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );


    //- Destructor
    virtual ~MultiComponentPhaseModel();


    // Member Functions

        //- Correct the phase properties other than the thermo
        void updateThermo();

        //- Correct the phase properties other than the thermo
        virtual void correct();

        //- Correct the species fractions
        virtual void correctSpecies();

        // Species

            //- Return whether the phase is pure (i.e., not multi-component)
            virtual bool pure() const;

            //- Return the species names
            virtual const wordList& species() const;

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const;

            //- Return the species densities
            virtual const PtrList<volScalarField>& Yrho() const;

            //- Return a species mass fraction by name
            virtual const volScalarField& Y(const word& name) const;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& YRef();


        // Thermophysical transport

            //- Return the density field
            virtual tmp<volScalarField> rho() const;

            //- Return the compressibility field
            virtual tmp<volScalarField> psi() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MultiComponentPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
