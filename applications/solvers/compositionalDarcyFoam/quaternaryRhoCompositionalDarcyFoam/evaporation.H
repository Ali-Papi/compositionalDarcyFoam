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
            volScalarField& XNa = phase.YRef()[2];
            volScalarField& XCl = phase.YRef()[3];
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
    N_Nal = 1e3*(*XNa)*(*rhoW)/Mw;
    NCl_l = 1e3*(*XCl)*(*rhoW)/Mw;
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
        bx[cellI] = (NT[cellI]-NCO2[cellI])/(PH[cellI]-1) + (NT[cellI]-NH2O[cellI])/(PP[cellI]-1);
		cx[cellI] = NT[cellI]*(N_Nal[cellI]+NCl_l[cellI])/((PH[cellI]-1)*(PP[cellI]-1));
        x1[cellI] = (-bx[cellI] + Foam::sqrt((sqr(bx[cellI]) - 4*cx[cellI])))/2 ;
		x2[cellI] = (-bx[cellI] - Foam::sqrt((sqr(bx[cellI]) - 4*cx[cellI])))/2;

        if ((sqr(bx[cellI]) - 4*cx[cellI])<0)
		    {
			    Info << "under the radical is negative, NL cannot be a complex number " << endl;
		    }
		    else if ( x1[cellI] > 0 && x2[cellI] <= 0)
		    {
			    NLflash[cellI] = x1[cellI];
			    Info << "x1[cellI] > 0 " << endl;
			    Info << "x1 =  " << x1[cellI] << endl;
			    Info << "NT[cellI] =  " << NT[cellI] << endl;
		    }
		    else if (x2[cellI] > 0 &&  x1[cellI] <= 0)
		    {
			    NLflash[cellI] = x2[cellI];
			    Info << " x2[cellI] > 0 " << endl;
			    Info << "x2 =  " << x2[cellI] << endl;
			    Info << "NT[cellI] =  " << NT[cellI] << endl;
		    }
		    else if ( x1[cellI] > 0 && x2[cellI] > 0)
		    {
			    Info << "both NL values are positive and this is not physical, NL can have only 1 value " << endl;
			    NLflash[cellI] = x1[cellI]; //just to give it something
			    Info << "x1 =  " << x1[cellI] << endl;
			    Info << "x2 =  " << x2[cellI] << endl;
			    Info << "NT[cellI] =  " << NT[cellI] << endl;
		    }
		    else
		    {
			    Info << "both NL values are non-positive and this means N is full gas" << endl;
		    }
		    
		
        //Info << " NLflash[cellI] = " << NLflash[cellI] << endl;
        //Info << " NT[cellI] = " << NT[cellI] << endl;

		if (NLflash[cellI] < NT[cellI]) // || NLflash[cellI] < NLmax[cellI] ) // || sigma2[cellI] < 1
		// || sigma1[cellI] >1)
		{
            //XCO2->write();
            Info << " yH2O[cellI] transport before flash =  " << YH2O->internalField()[cellI] << endl;
            Info << " yH2O[cellI+1] transport before flash =  " << YH2O->internalField()[cellI+1] << endl;
            Info << " yH2O[cellI+2] transport before flash =  " << YH2O->internalField()[cellI+2] << endl;
            NL[cellI]  = NLflash[cellI];

            rhoG->internalField()[cellI] = (p[cellI] * Mg.value() * 1e-3) / (R.value() * T.value());
            alphaG->internalField()[cellI] = ((NT[cellI] - NL[cellI]) * Mg.value() * 1e-3) / rhoG->internalField()[cellI];
            alphaW->internalField()[cellI] = porosity[cellI] - alphaG->internalField()[cellI];
            rhoW->internalField()[cellI] = (NL[cellI] * Mw.value() * 1e-3) / alphaW->internalField()[cellI];

			YCO2->internalField()[cellI] = alphaG->internalField()[cellI]*(PP[cellI]-1)/(PP[cellI]-PH[cellI]); // NCO2/(PH*NT+(NT-NL)*(1-PH));
			XCO2->internalField()[cellI] = alphaW->internalField()[cellI]*YCO2->internalField()[cellI]*PH[cellI]/alphaG->internalField()[cellI];
			XNa->internalField()[cellI] = alphaW->internalField()[cellI]*N_Nal[cellI]/NL[cellI];
			XCl->internalField()[cellI] = alphaW->internalField()[cellI]*NCl_l[cellI]/NL[cellI];
			XH2O->internalField()[cellI]= alphaW->internalField()[cellI]- XCO2->internalField()[cellI] - XNa->internalField()[cellI] - XCl->internalField()[cellI];
			YH2O->internalField()[cellI]= alphaG->internalField()[cellI]- YCO2->internalField()[cellI];

            Info << " yCO2[cellI] =  " << YCO2->internalField()[cellI] << endl;
            Info << " xCO2[cellI] =  " << XCO2->internalField()[cellI] << endl;
            Info << " xCO2[cellI+1] =  " << XCO2->internalField()[cellI+1] << endl;
            Info << " xCO2[cellI+2] =  " << XCO2->internalField()[cellI+2] << endl;
            Info << " yH2O[cellI] =  " << YH2O->internalField()[cellI] << endl;
            Info << " yH2O[cellI+1] =  " << YH2O->internalField()[cellI+1] << endl;
            Info << " yH2O[cellI+2] =  " << YH2O->internalField()[cellI+2] << endl;

            //Info << " XCO2[cellI] =  " << XCO2[cellI] << endl;
            //Info << " XCO2 =  " << XCO2 << endl;
            //Info << " cellI =  " << cellI << endl;
			
			
		}
		
		//Sb[cellI] = 1e-3*(NL[cellI]*Mw.value())/(rhob.value()*eps[cellI]);
		
	    }
    }
    
    //Sg = 1-Sb;
    
    
    //Info << "Sb evaporation: " << " Min(Sb) = " << gMin(Sb) << " Max(Sb) = " << gMax(Sb) << endl;

}
