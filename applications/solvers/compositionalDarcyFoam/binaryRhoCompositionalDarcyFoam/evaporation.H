{
    volScalarField& p = fluid.p();
    const volScalarField& porosity = fluid.porosity(); // Porosity
    PH = p/H;
    PP = p/Psat;


forAll(phases, i)
{
    phaseModel& phase = phases[i];

    volScalarField& alpha = phase;

    if (!phase.pure() && !phase.stationary())
    {
        if (phase.name() == "gas")
        {
            volScalarField& alphaG = phase;
            tmp<volScalarField> tRhoG = phase.rho();
            volScalarField& rhoG = const_cast<volScalarField&>(tRhoG());
            //rhoG = new volScalarField(tRhoG()); 
            //volScalarField& rhoG = phase.rho();          // Density
            volScalarField& YCO2 = phase.YRef()[0];
            volScalarField& YH2O = phase.YRef()[1];
        }
        else if (phase.name() == "water")
        {
            volScalarField& alphaW = phase;
            tmp<volScalarField> tRhoW = phase.rho();
            volScalarField& rhoW = const_cast<volScalarField&>(tRhoW());
            //rhoW = new volScalarField(tRhoW());
            //volScalarField& rhoW = phase.rho();          // Density
            volScalarField& XCO2 = phase.YRef()[0];
            volScalarField& XH2O = phase.YRef()[1];
        }
    }
}
Info << " rhoG->internalField()[cellI] bfroe =  " << rhoG->internalField()[0] << endl;
Info << " rhoW->internalField()[cellI] bfroe =  " << rhoW->internalField()[0] << endl;

    NH2O = 1e3*(*XH2O)*(*rhoW)/Mw+1e3*(*YH2O)*(*rhoG)/Mg;
    NCO2 = 1e3*(*XCO2)*(*rhoW)/Mw+1e3*(*YCO2)*(*rhoG)/Mg;
    NT = 1e3*(*alphaW)*(*rhoW)/Mw+1e3*(*alphaG)*(*rhoG)/Mg;
    //NCO2g = 1e3*eps*yCO2*rhoa*(1-Sb)/(Mg);
    NCO2l = 1e3*(*XCO2)*(*rhoW)/Mw;
    //NH2Ol = 1e3*eps*xH2O*rhob*Sb/(Mw);
    NH2Og = 1e3*(*YH2O)*(*rhoG)/Mg;
    NL = 1e3*(*alphaW)*(*rhoW)/(Mw);
    NCO2g = 1e3*(*YCO2)*(*rhoG)/Mg; //NCO2g = yCO2*(NT-NL);
    NH2Ol = 1e3*(*XH2O)*(*rhoW)/Mw; //NH2Ol = xH2O*NL;
    NG = NT-NL;
    //volScalarField sigma1 = NCO2/(NT*PH) + NH2O/(NT*PP) ;
    //volScalarField sigma2 = NCO2*PH/NT + NH2O*PP/NT ;
    volScalarField NLCO2min = NL;
    volScalarField NCO2_t = NCO2l+ NCO2g;
    NCO2g_t = NCO2g;
    NH2Ol_t = NH2Ol;


    // -------------------------------- method 2 ---------------------------   
    
    forAll(NCO2g, cellI)
    {
	    if (NCO2g[cellI]> 0)
	    {
		NLflash[cellI] = (NH2O[cellI])/(1-PH[cellI]) - (NCO2[cellI])/(PP[cellI]-1);
		
        Info << " NLflash[cellI] = " << NLflash[cellI] << endl;
        Info << " NT[cellI] = " << NT[cellI] << endl;

		if (NLflash[cellI] < NT[cellI]) // || NLflash[cellI] < NLmax[cellI] ) // || sigma2[cellI] < 1
		// || sigma1[cellI] >1)
		{
            //XCO2->write();
            Info << " alphaW->internalField()[cellI] before =  " << alphaW->internalField()[cellI] << endl;
            Info << " alphaG->internalField()[cellI] bfroe =  " << alphaG->internalField()[cellI] << endl;
            Info << " rhoG->internalField()[cellI] bfroe =  " << rhoG->internalField()[cellI] << endl;
            Info << " rhoW->internalField()[cellI] bfroe =  " << rhoW->internalField()[cellI] << endl;
			NL[cellI]  = NLflash[cellI];
            //dimensionedScalar unitPressure("unitPressure", dimPressure, 1.0);
            rhoG->internalField()[cellI] = (p[cellI] * Mg.value() * 1e-3) / (R.value() * T.value());
            //const volScalarField& rhoG = trhoG();
            //rhoG = (p*Mg*1e-3)/(R*T);
            //dimensionedScalar unitDensity("unitDensity", dimDensity, 1.0);
            //dimensionedScalar unitNumber("unitNumber", dimMoles/dimVolume, 1.0);  // [mol/m³]
            alphaG->internalField()[cellI] = ((NT[cellI] - NL[cellI]) * Mg.value() * 1e-3) / rhoG->internalField()[cellI];
            //const volScalarField& alphaG = tAlphaG();
            //alphaG = ((NT-NL)*Mg*1e-3)/rhoG;
            alphaW->internalField()[cellI] = porosity[cellI] - alphaG->internalField()[cellI];
            //const volScalarField& alphaW = tAlphaW();
            //alphaW = porosity - alphaG;
            //Sg = 1-Sb;

            rhoW->internalField()[cellI] = (NL[cellI] * Mw.value() * 1e-3) / alphaW->internalField()[cellI];
            //const volScalarField& rhoW = trhoW();
            Info << " alphaW->internalField()[cellI] after =  " << alphaW->internalField()[cellI] << endl;
            Info << " alphaG->internalField()[cellI] after =  " << alphaG->internalField()[cellI] << endl;
            Info << " rhoG->internalField()[cellI] after =  " << rhoG->internalField()[cellI] << endl;
            Info << " rhoW->internalField()[cellI] after =  " << rhoW->internalField()[cellI] << endl;

			YCO2->internalField()[cellI] = alphaG->internalField()[cellI]*(PP[cellI]-1)/(PP[cellI]-PH[cellI]);
			XCO2->internalField()[cellI] = alphaW->internalField()[cellI]*YCO2->internalField()[cellI]*PH[cellI]/alphaG->internalField()[cellI];
            Info << " yCO2[cellI] =  " << YCO2->internalField()[cellI] << endl;
            Info << " xCO2[cellI] =  " << XCO2->internalField()[cellI] << endl;
            Info << " xCO2[cellI+1] =  " << XCO2->internalField()[cellI+1] << endl;
            Info << " xCO2[cellI+2] =  " << XCO2->internalField()[cellI+2] << endl;

            //Info << " XCO2[cellI] =  " << XCO2[cellI] << endl;
            //Info << " XCO2 =  " << XCO2 << endl;
            //Info << " cellI =  " << cellI << endl;
			YH2O->internalField()[cellI]= alphaG->internalField()[cellI]- YCO2->internalField()[cellI];
			XH2O->internalField()[cellI]= alphaW->internalField()[cellI]- XCO2->internalField()[cellI];
		}
		
		//Sb[cellI] = 1e-3*(NL[cellI]*Mw.value())/(rhob.value()*eps[cellI]);
		
	    }
    }
    
    //Info << "Sb evaporation: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;

}
