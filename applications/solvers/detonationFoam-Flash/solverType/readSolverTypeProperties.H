word solverType(solverTypeProperties.lookupOrDefault<word>("solverType", "Euler"));

bool viscous(false);
word transportModelName("none");

if (solverType == "Euler")
{
    viscous = false;
    transportModelName = "none";
}
else if (solverType == "NS_Sutherland")
{
    viscous = true;
    transportModelName = "Sutherland";
}
else if (solverType == "NS_mixtureAverage")
{
    viscous = true;
    transportModelName = "mixtureAverage";
}
else if
(
    solverTypeProperties.found("viscous")
 || solverTypeProperties.found("transportModel")
)
{
    Info<< "Non-standard solverType, using explicit viscous/transportModel"
        << " entries" << endl;
}
else
{
    FatalIOErrorIn(args.executable().c_str(), solverTypeProperties)
        << "Unknown solverType " << solverType << nl
        << "Valid solverType entries:" << nl
        << "    Euler" << nl
        << "    NS_Sutherland" << nl
        << "    NS_mixtureAverage" << nl
        << "or set viscous/transportModel entries directly"
        << exit(FatalIOError);
}

viscous = Switch
(
    solverTypeProperties.lookupOrDefault<Switch>("viscous", Switch(viscous))
);
transportModelName = solverTypeProperties.lookupOrDefault<word>
(
    "transportModel",
    transportModelName
);

const bool mixtureAverageTransport(transportModelName == "mixtureAverage");
const bool sutherlandTransport(transportModelName == "Sutherland");

if
(
    transportModelName != "none"
 && transportModelName != "Sutherland"
 && transportModelName != "mixtureAverage"
)
{
    FatalIOErrorIn(args.executable().c_str(), solverTypeProperties)
        << "Unknown transportModel " << transportModelName << nl
        << "Valid transportModel entries: none, Sutherland, mixtureAverage"
        << exit(FatalIOError);
}

const Switch dropletsOn
(
    solverTypeProperties.lookupOrDefault<Switch>("droplets", Switch(false))
);
const word cloudName
(
    solverTypeProperties.lookupOrDefault<word>("cloudName", "sprayCloud")
);

Info<< "solverType = " << solverType << endl;
Info<< "viscous = " << Switch(viscous) << endl;
Info<< "transportModel = " << transportModelName << endl;
Info<< "droplets = " << dropletsOn << endl;
