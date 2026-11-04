/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    detonationFoam-Flash

Description
    Unified density-based compressible reacting flow solver for detonation
    simulations, based on central-upwind schemes of Kurganov and Tadmor with
    support for mesh-motion and topology changes.

    The solver merges the three former solver types (Euler, NS_Sutherland,
    NS_mixtureAverage) into a single executable.  The solver functionality is
    switched through the constant/solverTypeProperties dictionary:

        solverType       Euler | NS_Sutherland | NS_mixtureAverage (legacy)
        viscous          bool    - include viscous terms
        transportModel   none | Sutherland | mixtureAverage
        droplets         bool    - Lagrangian droplet (spray) cloud on/off

    Optional Lagrangian droplet cloud (basicSprayCloud) provides two-way
    coupling with the gas phase: evaporation mass source (continuity and
    species), droplet-gas drag and vapor momentum source (momentum), and
    convective + latent heat source (energy).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "fluidThermoMomentumTransportModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "psiReactionThermo.H"
#include "multivariateScheme.H"
#include "CombustionModel.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "fluxScheme.H"
#include "SLGThermo.H"
#include "basicSprayCloud.H"


int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"

    doubleScalar SW_position_temp;
    doubleScalar SW_position = 0;
    doubleScalar SW_position_limit = 0;

    turbulence->validate();
    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    scalar CoNum = 0.0;

    bool SW_track(solverType != "NS_mixtureAverage");
    SW_track = Switch(solverTypeProperties.lookupOrDefault<Switch>("SW_track", Switch(SW_track)));
    SW_position_limit = solverTypeProperties.lookupOrDefault<doubleScalar>
    (
        "SW_position_limit",
        1.0
    );
    Info<< "SW_track = " << Switch(SW_track) << endl;
    Info<< "SW_position_limit = "<< SW_position_limit << endl;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        if (mixtureAverageTransport)
        {
            #include "solverType/updateTransProperties.H"
        }

        #include "detoCellular.H"

        runTime++;

        if (dropletsOn)
        {
            parcelsPtr().storeGlobalPositions();
        }

        mesh.update();

        #include "calculateDt.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (dropletsOn)
        {
            parcelsPtr().evolve();
        }

        fluxSchemeFields->update(rho,rhoU,e,rPsi,c,phi,rhoPhi,rhoUPhi,rhoEPhi);

        #include "solverType/rhoEqn.H"
        #include "solverType/rhoUEqn.H"
        #include "solverType/rhoYEqn.H"
        #include "solverType/rhoEEqn.H"

        if (viscous)
        {
            turbulence->correct();
        }
        thermophysicalTransport->correct();

        runTime.write();

        if (SW_track)
        {
            for (label cellI=0; cellI<mesh.C().size(); cellI++)
            {
                if (p[cellI] > 101325+100)
                {
                    SW_position_temp = mesh.C()[cellI].x();
                    if (SW_position_temp > SW_position)
                    {
                        SW_position = SW_position_temp;
                    }
                }
            }
            reduce(SW_position,maxOp<scalar>());
            Info<< "Leading shock location = "<< SW_position << " m" << nl << endl;

            if (SW_position >= SW_position_limit) exit(EXIT_SUCCESS);
        }

        Info<< "ExecutionTime = "<< runTime.elapsedCpuTime() << " s" << nl << endl;
        Info<< "ClockTime = " << runTime.elapsedClockTime() << " s"<< nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

