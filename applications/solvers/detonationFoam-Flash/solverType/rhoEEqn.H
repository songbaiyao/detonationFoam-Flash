if (viscous && sutherlandTransport)
{
    solve
    (
        fvm::ddt(rhoE)
      + fvc::div(rhoEPhi)
      - fvc::div(sigmaDotU)
      ==
        reaction->Qdot()
    );
}
else
{
    solve
    (
        fvm::ddt(rhoE)
      + fvc::div(rhoEPhi)
      ==
        reaction->Qdot()
    );
}

e = rhoE/rho - 0.5*magSqr(U);
e.correctBoundaryConditions();
thermo.correct();

bool eCorrected = false;

if (mixtureAverageTransport)
{
    rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField()+0.5*magSqr(U.boundaryField()));

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, e) - fvc::ddt(rho, e)
        - fvm::laplacian(mixAlpha+turbulence->mut()/Prt, e)
        + diffsh
        - fvc::div(vk)
      ==
        reaction->Qdot()
    );

    if (dropletsOn)
    {
        EEqn -= parcelsPtr().Sh(e);
    }

    EEqn.relax();

    EEqn.solve("e");

    eCorrected = true;
}
else if (viscous)
{
    rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField()+0.5*magSqr(U.boundaryField()));

    if (dropletsOn)
    {
        solve
        (
            fvm::ddt(rho, e) - fvc::ddt(rho, e)
          + thermophysicalTransport->divq(e)
         ==
            parcelsPtr().Sh(e)
        );
    }
    else
    {
        solve
        (
            fvm::ddt(rho, e) - fvc::ddt(rho, e)
          + thermophysicalTransport->divq(e)
        );
    }

    eCorrected = true;
}
else if (dropletsOn)
{
    solve
    (
        fvm::ddt(rho, e) - fvc::ddt(rho, e)
     ==
        parcelsPtr().Sh(e)
    );

    eCorrected = true;
}

if (eCorrected)
{
    thermo.correct();
}

gamma=thermo.gamma();
R_gas=thermo.Cp()-thermo.Cv();

if (eCorrected)
{
    rhoE = rho*(e + 0.5*magSqr(U));
}
else
{
    rhoE.boundaryFieldRef() == rho.boundaryField()*(e.boundaryField()+0.5*magSqr(U.boundaryField()));
}

p.ref() = rho()/psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

Info<< "min/max(p) = "<< min(p).value() << ", " << max(p).value() << endl;
Info<< "min/max(T) = "<< min(T).value() << ", " << max(T).value() << endl;
