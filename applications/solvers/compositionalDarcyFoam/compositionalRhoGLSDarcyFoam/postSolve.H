// Retrieve solution
Ai = 0;
pAsBlock->retrieveSolution(Ai, p.internalField());

p.correctBoundaryConditions();

forAll(phases, solvePhasei)
{
    phaseModel& phase = phases[solvePhasei];

    if (phase.stationary()) continue;

    volScalarField& alpha = phase;

    if (phase.pure())
    {
        Ai++;

        pAsBlock->retrieveSolution(Ai, alpha.internalField());
        alpha.correctBoundaryConditions();
    }
    else
    {
        forAll(phase.Y(), compi)
        {
            Ai++;
            volScalarField& Yi = phase.YRef()[compi];

            pAsBlock->retrieveSolution(Ai, Yi.internalField());
            Yi.max(0);
            Yi.correctBoundaryConditions();

            if (compi == 0)
            {
                alpha = Yi;
            }
            else
            {
                alpha += Yi;
            }
        }
    }

    phase.alphafRef() = fvc::interpolate(alpha);
}

// Report the phase fractions and the phase fraction sum
// Calculate the void fraction
volScalarField alphaVoid
(
    IOobject
    (
        "alphaVoid",
        mesh.time().timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("one", dimless, 1)
);

// Calculate the porosity fraction
volScalarField porosity
(
    IOobject
    (
        "porosity",
        mesh.time().timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("one", dimless, 1)
);

forAll(fluid.stationaryPhases(), phasei)
{
    phaseModel& phase = fluid.stationaryPhases()[phasei];
    alphaVoid -= phase;
    if (phase.porous())
    {
        porosity -= phase;
    }
    
    Info<< phase.name() << " fraction, min, max = "
        << phase.weightedAverage(mesh.V()).value()
        << ' ' << min(phase).value()
        << ' ' << max(phase).value()
        << endl;
}

volScalarField sumAlphaMoving
(
    IOobject
    (
        "sumAlphaMoving",
        mesh.time().timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0)
);

forAll(fluid.movingPhases(), phasei)
{
    phaseModel& phase = fluid.movingPhases()[phasei];
    volScalarField& alpha = phase;

    Info<< phase.name() << " fraction, min, max = "
        << phase.weightedAverage(mesh.V()).value()
        << ' ' << min(phase).value()
        << ' ' << max(phase).value()
        << endl;
    
    alpha.max(0);
    sumAlphaMoving += alpha;
}

Info<< "Phase-sum volume fraction, min, max = "
    << (sumAlphaMoving + 1 - alphaVoid)()
      .weightedAverage(mesh.V()).value()
    << ' ' << min(sumAlphaMoving + 1 - alphaVoid).value()
    << ' ' << max(sumAlphaMoving + 1 - alphaVoid).value()
    << endl;

// Correct the sum of the phase fractions to avoid drift
forAll(fluid.movingPhases(), movingPhasei)
{
    fluid.movingPhases()[movingPhasei] *= alphaVoid/sumAlphaMoving;
}

// Update alphaVoid
fluid.alphaVoid() = alphaVoid;

// Update porosity
fluid.porosity() = porosity;

// Update Saturations
forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];

    if (phase.porous()) continue;

    const volScalarField& alphai = phase;
    phase.SRef() = alphai/fluid.porosity();
}

// snGrad(p_rgh)
fluid.snGradPSf(fvc::snGrad(p)*mesh.magSf());

// Correct fluxes and velocities
{
    surfaceScalarField alphaPhis
    (
        IOobject
        (
            "alphaPhis",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0), 0)
    );

    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];

        if (phase.pure())
        {
            label Ai = phase.blockIndex()[0];

            phase.alphaPhiRef() =
                eqnBlock[Ai][Ai].flux() + eqnBlock[Ai][0].flux() - fluid.phi0s()[Ai];
        }
        else
        {
            forAll(phase.Y(), compi)
            {
                label Ai = phase.blockIndex()[compi];

                if (compi == 0)
                {
                    phase.alphaPhiRef() =
                        eqnBlock[Ai][Ai].flux() + eqnBlock[Ai][0].flux() - fluid.phi0s()[Ai];
                }
                else
                {
                    phase.alphaPhiRef() +=
                        eqnBlock[Ai][Ai].flux() + eqnBlock[Ai][0].flux() - fluid.phi0s()[Ai];
                }
            }
        }

        phase.alphaRhoPhiRef() = 
            fvc::interpolate(phase.rho())*phase.alphaPhi();

        phase.phiRef() =
            phase.alphaPhi()/max(phase.alphaf(), phase.residualAlpha());

        alphaPhis += phase.alphaPhi();
    }

    phip = - alphaTaus*fluid.snGradPSf();
    phig = phiHbyA;
    phi = alphaPhis;

    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];

        phase.URef() = fvc::reconstruct(phase.alphaPhi())/max(phase, phase.residualAlpha());

        phase.URef().correctBoundaryConditions();
    }

    // Update fields and models if necessary
    fluid.update();
}

// updated densities
rho = fluid.rho();

