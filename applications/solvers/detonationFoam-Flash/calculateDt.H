rPsi = 1.0/psi;
c = sqrt(thermo.Cp()/thermo.Cv()*rPsi);

if (mesh.nInternalFaces())
{
    scalarField sumPhi(fvc::surfaceSum(mag(phi))().primitiveField());
    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
}
Info<< "Max Courant Number = " << " " << CoNum << endl;








