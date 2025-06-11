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
    Foam::phaseSystem

Description

SourceFiles
    phaseSystem.C
    newPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSystem_H
#define phaseSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "IOdictionary.H"
#include "phaseModel.H"
#include "phaseInterface.H"
#include "phaseInterfaceKey.H"
#include "convectionScheme.H"
#include "hashedWordList.H"
#include "HashPtrTable.H"
#include "runTimeSelectionTables.H"
#include "fixedValueFvsPatchFields.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseSystem
:
    public IOdictionary
{
public:

    // Public Typedefs

        typedef PtrList<fvScalarMatrix> eqnTable;

        typedef PtrList<eqnTable> eqnBlockTable;

        typedef PtrList<phaseModel> phaseModelList;

        typedef UPtrList<phaseModel> phaseModelPartialList;


protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Pressure field
        volScalarField p_;

        //- alphaVoid field
        volScalarField alphaVoid_;

        //- Total porosity field
        volScalarField porosity_;

        //- Flux snGrad(p)
        tmp<surfaceScalarField> snGradPSf_;

        //- Name of optional reference phase which is not solved for
        // but obtained from the sum of the other phases
        word referencePhaseName_;

        //- Phase models
        phaseModelList phaseModels_;

        //- Total volumetric flux
        surfaceScalarField phi_;

        //- Partial volumetric flux - Pressure
        surfaceScalarField phip_;

        //- Partial volumetric flux - Gravity
        surfaceScalarField phig_;

        //- Convection scheme for pressure linearisation
        tmp<fv::convectionScheme<scalar>> mvConvP_;

        //- Convection scheme for gravity contribution
        tmp<fv::convectionScheme<scalar>> mvConvG_;

        //- Flag for capillarity
        bool hasCapillarity_;

        //- Moving phase models
        phaseModelPartialList movingPhaseModels_;

        //- Stationary phase models
        phaseModelPartialList stationaryPhaseModels_;

        //- Porous phase models
        phaseModelPartialList porousPhaseModels_;

        //- Anisothermal phase models
        phaseModelPartialList anisothermalPhaseModels_;

        //- Multi-component phase models
        phaseModelPartialList multiComponentPhaseModels_;


    // Protected member functions

        //- Calculate and return the mixture flux
        tmp<surfaceScalarField> phiCalc() const;


public:

    //- Runtime type information
    TypeName("phaseSystem");

    //- Default name of the phase properties dictionary
    static const word propertiesName;


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseSystem,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from fvMesh
        phaseSystem(const fvMesh& mesh);


    // Selectors

        static autoPtr<phaseSystem> New
        (
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~phaseSystem();

    // Member Functions

        // Access

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return the pressure field
            inline const volScalarField& p() const;

            //- Access the pressure field
            inline volScalarField& p();

            //- Return the alphaVoid field
            inline const volScalarField& alphaVoid() const;

            //- Access the alphaVoid field
            inline volScalarField& alphaVoid();

            //- Return the porosity field
            inline const volScalarField& porosity() const;

            //- Access the porosity field
            inline volScalarField& porosity();

            //- Return the snGrad(p) (flux)
            inline tmp<surfaceScalarField> snGradPSf() const;

            //- Access the snGrad(p) (flux)
            inline void snGradPSf(tmp<surfaceScalarField> snGradPSf);

            //- Return the phase models
            inline const phaseModelList& phases() const;

            //- Access the phase models
            inline phaseModelList& phases();

            //- Return the models for phases that are moving
            inline const phaseModelPartialList& movingPhases() const;

            //- Access the models for phases that are moving
            inline phaseModelPartialList& movingPhases();

            //- Return the models for phases that are stationary
            inline const phaseModelPartialList& stationaryPhases() const;

            //- Access the models for phases that are stationary
            inline phaseModelPartialList& stationaryPhases();

            //- Return the models for phases that are porous
            inline const phaseModelPartialList& porousPhases() const;

            //- Access the models for phases that are porous
            inline phaseModelPartialList& porousPhases();

            //- Return the models for phases that have variable temperature
            inline const phaseModelPartialList& anisothermalPhases() const;

            //- Access the models for phases that have variable temperature
            inline phaseModelPartialList& anisothermalPhases();

            //- Return the models for phases that have multiple species
            inline const phaseModelPartialList& multiComponentPhases() const;

            //- Access the models for phases that have multiple species
            inline phaseModelPartialList& multiComponentPhases();

            //- Return the mixture flux
            inline const surfaceScalarField& phi() const;

            //- Access the mixture flux
            inline surfaceScalarField& phi();

            //- Return the partial flux - pressure
            inline const surfaceScalarField& phip() const;

            //- Access the partial flux - pressure
            inline surfaceScalarField& phip();

            //- Return the partial flux - gravity
            inline const surfaceScalarField& phig() const;

            //- Access the partial flux - gravity
            inline surfaceScalarField& phig();

            //- Return the convection scheme
            inline tmp<fv::convectionScheme<scalar>> mvConvP() const;

            //- Access the convection scheme
            inline void mvConvP(tmp<fv::convectionScheme<scalar>> mvConv);

            //- Return the convection scheme
            inline tmp<fv::convectionScheme<scalar>> mvConvG() const;

            //- Access the convection scheme
            inline void mvConvG(tmp<fv::convectionScheme<scalar>> mvConv);


        // PhaseList utilities

            //- Check availability of a phase
            bool foundPhase(const word& name) const;

            //- Return a phaseModel with a name
            const phaseModel& lookupPhase
            (
                const word& name
            ) const;


        // Sub-model construction

            //- Generate interfacial-model tables
            template<class ModelType>
            void generateInterfacialModels
            (
                const dictionary& dict,
                HashTable
                <
                    autoPtr<ModelType>,
                    phaseInterfaceKey,
                    phaseInterfaceKey::hash
                >& models
            ) const;


        // Properties

            //- Return the mixture density
            tmp<volScalarField> rho() const;

            //- Return the mixture velocity
            tmp<volVectorField> U() const;


        // Evolution

            //- Correct the fluid properties and models before solving equations
            virtual void correct();

            //- Update the fluid properties and models after solving equations
            virtual void update();

            //- Correct fixed-flux BCs to be consistent with the velocity BCs
            void correctBoundaryFlux();


        // Transfers

            //- Check availability of capillarity
            inline bool hasCapillarity() const;

            //- Set availability of capillarity
            inline void setCapillarity(bool capillarity);

            //- Return the Taus (rAUfs)
            virtual const surfaceScalarField& taus
            (
                const word& fieldName
            ) const = 0;

            //- Return the phiAlpha0P0s
            virtual PtrList<surfaceScalarField>& phi0s() = 0;

            //- Add the matrices for the coupled algorithm
            virtual void addTerms
            (
                phaseSystem::eqnBlockTable& eqnBlock
            ) = 0;


        // IO

            //- Read base phaseProperties dictionary
            virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
    #include "phaseSystemTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
