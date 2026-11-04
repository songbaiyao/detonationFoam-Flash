if (dropletsOn)
{
    solve
    (
        fvm::ddt(rho)
      + fvc::div(rhoPhi)
     ==
        parcelsPtr().Srho(rho)
    );
}
else
{
    solve(fvm::ddt(rho) + fvc::div(rhoPhi));
}
