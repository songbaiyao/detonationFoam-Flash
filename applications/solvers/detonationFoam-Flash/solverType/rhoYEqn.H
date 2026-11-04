tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(rhoPhi,Yi_h)")
    )
);

if (mixtureAverageTransport)
{
    vk *= 0;
    diffsh *= 0;
    DiffError *= 0;
    surfaceScalarField sumOfOtherSpeices = fvc::interpolate(DiffError) & mesh.Sf();

    forAll(Y,i)
    {
        volScalarField& Yi = Y[i];
        volScalarField& Xi = X[i];
        volScalarField& Diffi = Diff[i];
        DiffError += (rho*Diffi*fvc::grad(Yi) + Yi*rho*Diffi*fvc::grad(MW)/MW);
        if(thermalDiffusion)
        {
            if(Y[i].name() == "H")
            {
               DiffError += (-rho*Diffi*TDRatio_H*fvc::grad(T)*Yi/(T*(Xi+SMALL)));
            }
            else if(Y[i].name() == "H2")
            {
               DiffError +=  (-rho*Diffi*TDRatio_H2*fvc::grad(T)*Yi/(T*(Xi+SMALL)));
            }
        }
    }

    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);
    surfaceScalarField  phiUc = fvc::interpolate(DiffError) & mesh.Sf();

    for(label i=0; i < Y.size(); i++)
    {
        if(i != inertIndex)
        {
            volScalarField& Yi = Y[i];
            volScalarField& Xi = X[i];
            volScalarField& Diffi = Diff[i];
            volScalarField  ThermalDiffRatio = TDRatio_H*0;

            if(thermalDiffusion)
            {
                if(Y[i].name() == "H")
                {
                   ThermalDiffRatio = TDRatio_H;
                }
                else if(Y[i].name() == "H2")
                {
                   ThermalDiffRatio = TDRatio_H2;
                }
                else
                {
                   ThermalDiffRatio *= 0;
                }
            }

            Hsi = composition.Hs(i, p, T);
            HEi = composition.HE(i, p, T);

            surfaceScalarField speicesFlux = fvc::interpolate(rho*Diffi*fvc::grad(Yi)+Yi*rho*Diffi*fvc::grad(MW)/MW+(-rho*Diff[i]*ThermalDiffRatio*fvc::grad(T)*Yi/(T*(Xi+SMALL)))-Yi*DiffError)& mesh.Sf();
            vk += fvc::interpolate(Hsi)*speicesFlux;
            sumOfOtherSpeices += speicesFlux;
            diffsh += fvc::laplacian(mixAlpha*HEi, Yi);
            Ydiffuision = rho*Diffi*Yi/MW*fvc::grad(MW);

            fvScalarMatrix YEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(rhoPhi, Yi)
              - fvm::laplacian(rho*Diffi, Yi)
              + mvConvection->fvmDiv(phiUc,Yi)
              - fvc::div(Ydiffuision)
              - fvc::laplacian(-rho*Diffi*ThermalDiffRatio*Yi/(T*(Xi+SMALL)),T)
              ==
                reaction->R(Yi)
            );

            if (dropletsOn)
            {
                YEqn -= parcelsPtr().SYi(i, Yi);
            }

            YEqn.relax();

            YEqn.solve("Yi");

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    {
        label i = inertIndex;

        Hsi = composition.Hs(i, p, T);
        HEi = composition.HE(i, p, T);

        vk += fvc::interpolate(Hsi)*-sumOfOtherSpeices;
        diffsh += fvc::laplacian(mixAlpha*HEi, 1-Yt);
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
else
{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(rhoPhi, Yi)
             ==
                reaction->R(Yi)
            );

            if (viscous)
            {
                YiEqn += thermophysicalTransport->divj(Yi);
            }

            if (dropletsOn)
            {
                YiEqn -= parcelsPtr().SYi(i, Yi);
            }

            YiEqn.relax();

            YiEqn.solve("Yi");

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
