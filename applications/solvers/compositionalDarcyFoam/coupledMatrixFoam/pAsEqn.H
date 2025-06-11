// Correct boundaries
forAll(phases, phasei)
{
    phaseModel& phase = phases[phasei];
    volScalarField& alpha = phase;

    alpha.correctBoundaryConditions();

    if (!phase.stationary())
    {
        phase.URef().correctBoundaryConditions();
    }

    if (!phase.pure())
    {
        forAll(phase.Y(), Yi)
            phase.YRef()[Yi].correctBoundaryConditions();
    }
}

// Correct fixed-flux BCs to be consistent with the velocity BCs
fluid.correctBoundaryFlux();

word phipScheme("div(phip,alpha)");
word phigScheme("div(phig,alpha)");

fluid.mvConvP
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        alphaTable,
        phip,
        mesh.schemesDict().divScheme(phipScheme)
    )
);

fluid.mvConvG
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        alphaTable,
        phig,
        mesh.schemesDict().divScheme(phigScheme)
    )
);

tmp<fv::convectionScheme<scalar>> mvConvG(fluid.mvConvG());

// Update phase system models
fluid.correct();

// Create all matrices for equations and cross-coupling terms
autoPtr<phaseSystem::eqnBlockTable> eqnBlockPtr
(
    new phaseSystem::eqnBlockTable()
);

phaseSystem::eqnBlockTable& eqnBlock = eqnBlockPtr();

eqnBlock.setSize(blockSize);

forAll(eqnBlock, blocki)
{
    eqnBlock.set(blocki, new phaseSystem::eqnTable());

    eqnBlock[blocki].setSize(blockSize);
}

// Total predicted gravity flux
surfaceScalarField phiHbyA
(
    IOobject
    (
        "phiHbyA",
        mesh.time().timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0), 0)
);

// Gravity flux
PtrList<surfaceScalarField> phigs(phases.size());
{
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        if (phase.stationary()) continue;

        const word tauName(IOobject::groupName("alphaTaug", phase.name()));

        surfaceScalarField phigi(fluid.taus(tauName)*fvc::interpolate(phase.rho())*gSf);

        phigs.set
        (
            phasei,
            phigi
        );

        forAll(phase.URef().boundaryField(), patchi)
        {
            if (phase.URef().boundaryField()[patchi].fixesValue())
            {
                phigs[phasei].boundaryField()[patchi] == Zero;
            }
        }

        phiHbyA += phigs[phasei];
    }
}

const surfaceScalarField& alphaTaus = fluid.taus("alphaTausp");

// Update the fixedFluxPressure BCs to ensure flux consistency
{
    volScalarField::GeometricBoundaryField& pb = p.boundaryField();

    forAll(pb, patchi)
    {
        if (isA<fluxPressureFvPatchScalarField>(pb[patchi]))
        {
            scalarField phib(pb[patchi].size(), 0);

            forAll(phases, phasei)
            {
                phaseModel& phase = phases[phasei];

                if (phase.stationary()) continue;

                phib +=
                    phase.boundaryField()[patchi]*phase.phi()().boundaryField()[patchi];
            }

            scalarField snGrad
            (
                (phiHbyA.boundaryField()[patchi] - phib)
              / (mesh.magSf().boundaryField()[patchi]*alphaTaus.boundaryField()[patchi])
            );

            refCast<fluxPressureFvPatchScalarField>(pb[patchi]).updateCoeffs(snGrad);
        }
    }
}

// Update snGrad(p) in phaseSystem
fluid.snGradPSf(fvc::snGrad(p)*mesh.magSf());

// Construct the transport part of the pressure equation
eqnBlock[0].set
(
    0,
  - fvm::laplacian(alphaTaus, p)
);

fvScalarMatrix& pEqn = eqnBlock[0][0]; 

pEqn.setReference(pRefCell, pRefValue);

forAll(phases, solvePhasei)
{
    phaseModel& phase = phases[solvePhasei];

    if (phase.stationary()) continue;

    volScalarField& alpha = phase;

    if (phase.pure())
    {
        Info<< " This line is being read from the pure loop. Even though we have species, my loop still is operated" << endl;
        // phase flux by gravity
        surfaceScalarField phigi
        (
            phigs[phase.index()]
          / max(mvConvG->interpolate(phig, alpha), phase.residualAlpha())
        );

        label Ai = phase.blockIndex()[0];

        const volScalarField& rhoi = phase.rho()();

        fvScalarMatrix divPhiG
        (
            mvConvG->fvmDiv(phigi, alpha)
          + fvm::Sp(((fvc::ddt(rhoi) + (phase.U() & fvc::grad(rhoi)))/rhoi), alpha)
          - fvc::Sp(fvc::ddt(rhoi)/rhoi, alpha)
        );

        // alpha in alphaEqn [Ai, Ai]
        eqnBlock[Ai].set
        (
            Ai,
            fvm::ddt(alpha)
          + divPhiG
        );

        // alpha in pEqn [p, Ai]
        eqnBlock[0].set(Ai, divPhiG);

        fvScalarMatrix rhoInp
        (
            (alpha/rhoi)*fvm::ddt(phase.psi(), p)
        );

        // [p, p]
        pEqn += rhoInp;

        // [Ai, p]
        eqnBlock[Ai].set(0, rhoInp);
    }
    else
    {
        Info<< " This is species loop (else)" << endl;
        // phase flux by gravity
        surfaceScalarField phigi
        (
            phigs[phase.index()]
          / max(mvConvG->interpolate(phig, alpha), phase.residualAlpha())
        );

        forAll(phase.Y(), compi)
        {
            volScalarField& Yi = phase.YRef()[compi];

            label Ai = phase.blockIndex()[compi];

            fvScalarMatrix divPhiG
            (
                mvConvG->fvmDiv(phigi, Yi)
            );

            // alpha in alphaEqn [Ai, Ai]
            eqnBlock[Ai].set
            (
                Ai,
                fvm::ddt(Yi)
              + divPhiG
            );

            // Yi in pEqn [p, Yi]
            eqnBlock[0].set
            (
                Ai,
                divPhiG
            );
        }
/*        const word tauName(IOobject::groupName("alphaTaug", phase.name()));

        surfaceScalarField phigi
        (
            fluid.taus(tauName)*gSf
          / max(mvConvG->interpolate(phig, alpha), phase.residualAlpha())
        );

        forAll(phase.URef().boundaryField(), patchi)
        {
            if (phase.URef().boundaryField()[patchi].fixesValue())
            {
                phigi.boundaryField()[patchi] == Zero;
            }
        }

        const label nComps = phase.Y().size();

        forAll(phase.Y(), compi)
        {
            volScalarField& Yi = phase.YRef()[compi];

            label Ai = phase.blockIndex()[compi];

            fvScalarMatrix divPhiG
            (
                mvConvG->fvmDiv(phigi*fvc::interpolate(phase.Yrho()[compi]), Yi)
            );

            forAll(phase.Y(), compj)
            {
                label Aj = phase.blockIndex()[compj];

                if (Ai == Aj)
                {
                    // alpha in alphaEqn [Ai, Ai]
                    eqnBlock[Ai].set
                    (
                        Ai,
                        fvm::ddt(Yi)
                      + divPhiG
                    );

                    // Yi in pEqn [p, Yi]
                    eqnBlock[0].set
                    (
                        Ai,
                        nComps*divPhiG
                    );
                }
                else
                {
                    // alpha in alphaEqn [Aj, Ai]
                    eqnBlock[Aj].set
                    (
                        Ai,
                        divPhiG
                    );
                }
            }
        }*/
    }
}

fluid.addTerms(eqnBlock);

forAll(eqnBlock, varI)
{
    if (!eqnBlock.set(varI)) continue;

    forAll(eqnBlock[varI], varJ)
    {
        if (!eqnBlock[varI].set(varJ)) continue;

        if (varI == varJ)
        {
            pAsBlock->insertEquation(varI, eqnBlock[varI][varJ]);
        }
        else
        {
            pAsBlock->insertEquationCoupling(varI, varJ, eqnBlock[varI][varJ]);
        }
    }
}

