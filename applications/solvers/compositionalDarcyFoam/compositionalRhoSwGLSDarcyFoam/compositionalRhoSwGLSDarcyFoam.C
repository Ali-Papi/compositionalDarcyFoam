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

Application
    compositionalSwDarcyFoam

Description
    The code for salt precipitation in CO2 storage in saline aquifers with interphase interaction in compositional multiphase flow.
    The solver only evaluates Gas-Liquid-solid (GLS) equilibrium conditions.
    The water saturation in the equilibrium stage is updated during flash and is not assumed constant.
    Multiphase solver for porous media, with implicit coupling between 
    pressure, phase fractions and compositions achieved by fvBlockMatrix.

Authors
    Ali Papi

\*---------------------------------------------------------------------------*/

#include "phaseSystem.H"
#include "fvBlockMatrix.H"
#include "blockSystem.H"
#include "fluxPressureFvPatchScalarField.H"
#include "multivariateSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControls.H"
    #include "createFields.H"

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readBlockSolverControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
       
        // --- Pressure-alphas corrector loop
        while (pAsBlock->loop())
        {
            // Initialize the pAs block system (matrix, source and reference to pAs)
            pAsBlock->newMatrix();

            // Assemble and insert equations
            #include "pAsEqn.H"

            // Solve the block matrix
            pAsBlock->solve();

            // Retrieve solution and post-solve
            #include "postSolve.H"
        }

        #include "evaporation.H"

        runTime.write();
        if (runTime.outputTime())
        {
            Ns.write(); 
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
