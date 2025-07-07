{
    volScalarField& porosity = fluid.porosity(); // Porosity
 
    delNs = Ns - NsLast;
    delPhi = delNs*Ms*1e-3/Rhos;
    porosity = porosity + delPhi;
}
