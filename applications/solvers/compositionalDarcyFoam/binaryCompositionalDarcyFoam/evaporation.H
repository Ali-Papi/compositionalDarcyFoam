{
    volScalarField& p = fluid.p();
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
            rhoG = new volScalarField(tRhoG()); 
            //volScalarField& rhoG = phase.rho();          // Density
            volScalarField& YCO2 = phase.YRef()[0];
            volScalarField& YH2O = phase.YRef()[1];
        }
        else if (phase.name() == "water")
        {
            volScalarField& alphaW = phase;
            tmp<volScalarField> tRhoW = phase.rho();
            rhoW = new volScalarField(tRhoW());
            //volScalarField& rhoW = phase.rho();          // Density
            volScalarField& XCO2 = phase.YRef()[0];
            volScalarField& XH2O = phase.YRef()[1];
        }
    }
}

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
			NL[cellI]  = NLflash[cellI];
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
    
    //Sg = 1-Sb;
    
    
    //Info << "Sb evaporation: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;

}
