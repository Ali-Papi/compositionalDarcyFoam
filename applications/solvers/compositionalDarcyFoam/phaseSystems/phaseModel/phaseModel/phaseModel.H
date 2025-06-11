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
    Foam::phaseModel

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseSystem;

/*---------------------------------------------------------------------------*\
                         Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField
{
    // Private Data

        //- Reference to the phaseSystem to which this phase belongs
        const phaseSystem& fluid_;

        //- Name of phase
        word name_;

        //- Index of phase
        label index_;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilise the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha_;

        //- Optional maximum phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        //- Dilatation rate
        //tmp<volScalarField> divU_;

        //- Phase-fraction interpolated
        surfaceScalarField alphaf_;

        //- Allocation in blockSystem
        labelList blockIndex_;


public:

    //- Runtime type information
    ClassName("phaseModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseModel,
            phaseSystem,
            (
                const phaseSystem& fluid,
                const word& phaseName,
                const bool referencePhase,
                const label index
            ),
            (fluid, phaseName, referencePhase, index)
        );


    // Constructors

        phaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );

        //- Return clone
        autoPtr<phaseModel> clone() const;


    // Selectors

        static autoPtr<phaseModel> New
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const bool referencePhase,
            const label index
        );

        //- Return a pointer to a new phase created on freestore
        //  from Istream
        class iNew
        {
            const phaseSystem& fluid_;
            const word& referencePhaseName_;
            mutable label indexCounter_;

        public:

            iNew
            (
                const phaseSystem& fluid,
                const word& referencePhaseName
            )
            :
                fluid_(fluid),
                referencePhaseName_(referencePhaseName),
                indexCounter_(-1)
            {}

            autoPtr<phaseModel> operator()(Istream& is) const
            {
                ++indexCounter_;

                const word phaseName(is);

                return autoPtr<phaseModel>
                (
                    phaseModel::New
                    (
                        fluid_,
                        phaseName,
                        phaseName == referencePhaseName_,
                        indexCounter_
                    )
                );
            }
        };


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        //- Return the name of this phase
        const word& name() const;

        //- Return the name of the phase for use as the keyword in PtrDictionary
        const word& keyword() const;

        //- Return the index of the phase
        label index() const;

        //- Return the system to which this phase belongs
        const phaseSystem& fluid() const;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilise the phase momentum as the phase-fraction -> 0
        const dimensionedScalar& residualAlpha() const;

        //- Return the maximum phase-fraction (e.g. packing limit)
        scalar alphaMax() const;

        //- Acess the phase-fraction interpolated
        surfaceScalarField& alphafRef();

        //- Return the phase-fraction interpolated
        tmp<surfaceScalarField> alphaf() const;

        //- Correct the phase properties
        virtual void correct();

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Correct the reactions
        virtual void correctReactions();

        //- Correct the species concentrations
        virtual void correctSpecies();

        //- Return the index of the phase in blockSystem
        const labelList& blockIndex() const;

        //- Acess the index of the phase in blockSystem
        labelList& blockIndex();

        //- Ensure that the flux at inflow/outflow BCs is preserved
        void correctInflowOutflow(surfaceScalarField& alphaPhi) const;

        //- Read phase properties dictionary
        virtual bool read();


        // Density variation and compressibility

            //- Return true if the phase is incompressible otherwise false
            //virtual bool incompressible() const = 0;

            //- Return true if the phase is constant density otherwise false
            //virtual bool isochoric() const = 0;

            //- Return the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            tmp<volScalarField> divU() const;

            //- Set the phase dilatation rate (d(alpha)/dt + div(alpha*phi))
            void divU(tmp<volScalarField> divU);


        // Thermo

            //- Return the thermophysical model
            //virtual const basicRhoThermo& thermo() const = 0;

            //- Access the thermophysical model
            //virtual basicRhoThermo& thermoRef() = 0;

            //- Return the density field
            virtual tmp<volScalarField> rho() const = 0;

            //- Return the compressibility field
            virtual tmp<volScalarField> psi() const = 0;

            //- Return whether the phase is isothermal
            virtual bool isothermal() const = 0;

            //- Return the enthalpy equation
            //virtual tmp<fvScalarMatrix> heEqn() = 0;


        // Species

            //- Return whether the phase is pure (i.e., not multi-component)
            virtual bool pure() const = 0;

            //- Return the species names
            virtual const wordList& species() const = 0;

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const = 0;

            //- Return the species densities
            virtual const PtrList<volScalarField>& Yrho() const = 0;

            //- Return a species mass fraction by name
            virtual const volScalarField& Y(const word& name) const = 0;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& YRef() = 0;

            //- Return the fuel consumption rate matrix
            virtual tmp<fvScalarMatrix> R(volScalarField& Yi) const = 0;


        // Momentum

            //- Return whether the phase is stationary
            virtual bool stationary() const = 0;

            //- Return whether the phase is porous
            virtual bool porous() const = 0;

            //- Return the velocity
            virtual tmp<volVectorField> U() const = 0;

            //- Access the velocity
            virtual volVectorField& URef() = 0;

            //- Return the saturation
            virtual tmp<volScalarField> S() const = 0;

            //- Access the saturation
            virtual volScalarField& SRef() = 0;

            //- Return the absolut permeability
            virtual tmp<volDiagTensorField> invK() const = 0;

            //- Access the absolut permeability
            virtual volDiagTensorField& invKRef() = 0;

            //- Return the absolut permeability vectorField at faces
            virtual tmp<surfaceDiagTensorField> invKTf() const = 0;

            //- Return the absolut permeability at faces
            virtual tmp<surfaceScalarField> invKf() const = 0;

            //- Return the volumetric flux
            virtual tmp<surfaceScalarField> phi() const = 0;

            //- Access the volumetric flux
            virtual surfaceScalarField& phiRef() = 0;

            //- Return the volumetric flux of the phase
            virtual tmp<surfaceScalarField> alphaPhi() const = 0;

            //- Access the volumetric flux of the phase
            virtual surfaceScalarField& alphaPhiRef() = 0;

            //- Return the mass flux of the phase
            virtual tmp<surfaceScalarField> alphaRhoPhi() const = 0;

            //- Access the mass flux of the phase
            virtual surfaceScalarField& alphaRhoPhiRef() = 0;

        // Transport

            //- Return the laminar dynamic viscosity
            virtual tmp<volScalarField> mu() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
