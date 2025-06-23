    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    uniformDimensionedScalarField hRef
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar ghRef(- mag(g)*hRef);
    volScalarField gh("gh", (g & mesh.C()) - ghRef);
    surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

    surfaceScalarField gSf(g & mesh.Sf());

    autoPtr<phaseSystem> fluidPtr
    (
        phaseSystem::New(mesh)
    );

    phaseSystem& fluid = fluidPtr();
    phaseSystem::phaseModelList& phases = fluid.phases();

    volScalarField& p = fluid.p();

    mesh.schemesDict().setFluxRequired(p.name());

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh
        ),
        fluid.rho()
    );

    // Create fields ref
    surfaceScalarField& phi = fluid.phi();
    surfaceScalarField& phip = fluid.phip();
    surfaceScalarField& phig = fluid.phig();

    dictionary blockSolverDict = mesh.solutionDict().subDict("blockSolver");

    label Ai = 0;
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        if (!phase.stationary() && phase.pure())
        {
            Ai++;
            phase.blockIndex().append(Ai);
        }
        else if (!phase.stationary() && !phase.pure())
        {
            forAll(phase.Y(), compi)
            {
                Ai++;
                phase.blockIndex().append(Ai);
            }
        }
    }
    const label blockSize = ++Ai;

    // set size of phi0s
    fluid.phi0s().setSize(blockSize);

    // Block vector field for pressure (first entry) and phases fractions 
    // (others entries).   
    autoPtr<blockSystem> pAsBlock(blockSystem::New(blockSize, pAsDict, mesh)); 

    //- Phase-fraction field table for multivariate discretisation
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable alphaTable;

    // Initialise p in pAs
    pAsBlock->blockAdd(0, p.internalField());

    // Update fields
    // Face volume fractions
    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];
        const volScalarField& alpha = phase;

        phase.alphafRef() = fvc::interpolate(alpha);

        if (phase.stationary()) continue;

        if (phase.pure())
        {
            alphaTable.add(alpha);

            // Initialise alpha in pAs
            label Ai = phase.blockIndex()[0];
            pAsBlock->blockAdd(Ai, alpha.internalField());
        }
        else
        {
            forAll(phase.Y(), compi)
            {
                const volScalarField& Y = phase.Y()[compi];

                alphaTable.add(Y);

                // Initialise Yi in pAs
                label Ai = phase.blockIndex()[compi];
                pAsBlock->blockAdd(Ai, Y.internalField());
                mesh.schemesDict().setFluxRequired(Y.name());
            }
        }
    }


    // ---------------------------------------------------------
    // Read the top-level dictionary
    IOdictionary thermophysicalProperties
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dictionary subDict0 = thermophysicalProperties.subDict("solutionSpecies");
    dictionary subDict1 = subDict0.subDict("Species");
    //dimensionedScalar D1 = subDict1.lookup("D1");
    //dimensionedScalar D2 = subDict1.lookup("D2"); //both phases should have it, if any.
    dimensionedScalar H("H", subDict1.lookup("H"));
    dimensionedScalar Mw("Mw", subDict1.lookup("Mw"));
    dimensionedScalar Mg("Mg", subDict1.lookup("Mg"));
    dimensionedScalar Psat("Psat", subDict1.lookup("Psat"));
    
    volScalarField PH = p/H;
    volScalarField PP = p/Psat;
    
    //const volScalarField& porosity = fluid.porosity();// Porosity

    // Declare as pointers (just declared, not initialized)
    // Declare actual volScalarField references instead of pointers
    volScalarField *alphaG = nullptr, *alphaW = nullptr;
    volScalarField *rhoG = nullptr, *rhoW = nullptr;
    volScalarField *YCO2 = nullptr, *YH2O = nullptr;
    volScalarField *XCO2 = nullptr, *XH2O = nullptr;
    volScalarField *YNa = nullptr, *YCl = nullptr;
    volScalarField *XNa = nullptr, *XCl = nullptr;


forAll(phases, i)
{
    phaseModel& phase = phases[i];
    if (!phase.pure() && !phase.stationary())
    {
        if (phase.name() == "gas")
        {
            alphaG = const_cast<volScalarField*>(&static_cast<const volScalarField&>(phase));
            tmp<volScalarField> tRhoG = phase.rho();
            rhoG = new volScalarField(tRhoG()); // or const volScalarField& rhoG = tRhoG();
            YCO2 = &phase.YRef()[0];
            YH2O = &phase.YRef()[1];
        }
        else if (phase.name() == "water")
        {
            alphaW = const_cast<volScalarField*>(&static_cast<const volScalarField&>(phase));
            tmp<volScalarField> tRhoW = phase.rho();
            rhoW = new volScalarField(tRhoW());
            XCO2 = &phase.YRef()[0];
            XH2O = &phase.YRef()[1];
            XNa = &phase.YRef()[2];
            XCl = &phase.YRef()[3];
        }
    }
}

// Check all fields are properly assigned
if (!alphaG || !rhoG || !YCO2 || !YH2O || !alphaW || !rhoW || !XCO2 || !XH2O || !XNa || !XCl)
{
    FatalErrorInFunction
        << "Phase variables not initialized properly. Check phase names or phase definitions."
        << exit(FatalError);
}

// Now compute derived fields
volScalarField NH2O = 1e3*(*XH2O)*(*rhoW)/Mw + 1e3*(*YH2O)*(*rhoG)/Mg;
volScalarField NCO2 = 1e3*(*XCO2)*(*rhoW)/Mw + 1e3*(*YCO2)*(*rhoG)/Mg;
//volScalarField NT = 1e3*(*alphaW)*(*rhoW)/Mw + 1e3*(*alphaG)*(*rhoG)/Mg + Ns;
volScalarField NL = 1e3*(*alphaW)*(*rhoW)/Mw;
volScalarField NG = 1e3*(*alphaG)*(*rhoG)/Mg; // NT - NL - Ns;
Info << nl << "Reading moles of salt, Ns" << endl;
    volScalarField Ns
    (
        IOobject
        (
            "Ns",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << nl << "Reading equilibrium constant, keq" << endl;
    volScalarField keq
    (
        IOobject
        (
            "keq",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

volScalarField NT = NL + NG + Ns;

volScalarField NH2Ol = 1e3*(*XH2O)*(*rhoW)/Mw;
volScalarField NH2Og = 1e3*(*YH2O)*(*rhoG)/Mg;
volScalarField NCO2l = 1e3*(*XCO2)*(*rhoW)/Mw;
volScalarField NCO2g = 1e3*(*YCO2)*(*rhoG)/Mg;
volScalarField N_Nal = 1e3*(*XNa)*(*rhoW)/Mw;
volScalarField NCl_l = 1e3*(*XCl)*(*rhoW)/Mw;
volScalarField NCl = NCl_l + 0.5*Ns;
volScalarField N_Na = N_Nal + 0.5*Ns;
volScalarField NLflash = NL;
volScalarField Nsflash = Ns;
volScalarField bx = (NT-NCO2-Ns)/(PH-1) + (NT-NH2O-Ns)/(PP-1);
volScalarField cx = (NT-Ns)*(N_Na+NCl-Ns)/((PH-1)*(PP-1));
volScalarField x1 = (-bx + sqrt((sqr(bx) - 4*cx)))/2;
volScalarField x2 = (-bx - sqrt((sqr(bx) - 4*cx)))/2;
volScalarField Bx = -2*(N_Na+NCl);
volScalarField Cx = 4*N_Na*NCl - 4*keq*sqr(NLflash)*sqr(Mw)*1e-6;
volScalarField X1 = (-Bx + sqrt((sqr(Bx) - 4*Cx)))/2;
volScalarField X2 = (-Bx - sqrt((sqr(Bx) - 4*Cx)))/2;

volScalarField NCO2g_t = NCO2g;
volScalarField NH2Ol_t = NH2Ol;

double NSNL_error = 1;
double equil_flag = 1;
volScalarField NsflashLastItr = Nsflash;
