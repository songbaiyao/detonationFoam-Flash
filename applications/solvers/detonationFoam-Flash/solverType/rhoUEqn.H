if (viscous && !mixtureAverageTransport)
{
    muEff = turbulence->muEff();
    tauMC = muEff*dev2(Foam::T(fvc::grad(U)));
}

solve(fvm::ddt(rhoU) + fvc::div(rhoUPhi));

U.ref() = rhoU()/rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

if (mixtureAverageTransport)
{
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) - fvc::ddt(rho, U)
      - fvm::laplacian(mixMu+turbulence->mut(),U)
      - fvc::div((mixMu+turbulence->mut())*dev2(Foam::T(fvc::grad(U))))
    );

    if (dropletsOn)
    {
        UEqn -= parcelsPtr().SU(U);
    }

    UEqn.relax();

    UEqn.solve("U");
}
else if (viscous)
{
    if (dropletsOn)
    {
        solve
        (
            fvm::ddt(rho, U) - fvc::ddt(rho, U)
          - fvm::laplacian(muEff, U)
          - fvc::div(tauMC)
         ==
            parcelsPtr().SU(U)
        );
    }
    else
    {
        solve
        (
            fvm::ddt(rho, U) - fvc::ddt(rho, U)
          - fvm::laplacian(muEff, U)
          - fvc::div(tauMC)
        );
    }
}
else if (dropletsOn)
{
    solve
    (
        fvm::ddt(rho, U) - fvc::ddt(rho, U)
     ==
        parcelsPtr().SU(U)
    );
}

if (viscous && sutherlandTransport)
{
    sigmaDotU =
    (
        (
            fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
          + fvc::dotInterpolate(mesh.Sf(), tauMC)
        )
      & fvc::interpolate(U)
    );
}

if (viscous || dropletsOn)
{
    rhoU = rho*U;
}
