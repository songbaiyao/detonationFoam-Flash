Switch differentialDiffusion(false);
Switch thermalDiffusion(false);
dimensionedScalar Prt("Prt", dimless, 0.7);
dimensionedScalar Sct("Sct", dimless, 0.7);

if (mixtureAverageTransport)
{
    Info<< "Reading chemistry properties\n" << endl;

    IOdictionary chemistryProperties
    (
        IOobject
        (
            "chemistryProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    differentialDiffusion = Switch
    (
        chemistryProperties.lookupOrDefault("differentialDiffusion", true)
    );
    thermalDiffusion = Switch
    (
        chemistryProperties.lookupOrDefault("thermalDiffusion", false)
    );
    Prt = dimensionedScalar(chemistryProperties.lookup("Prt"));
    Sct = dimensionedScalar(chemistryProperties.lookup("Sct"));
}


volScalarField MW
(
    IOobject
    (
        "MW",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("MW", dimensionSet (0,0,0,0,0,0,0),0.0)
);

volVectorField Ydiffuision
(
    IOobject
    (
        "Ydiffuision",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("Ydiffuision", dimDensity*dimVelocity*dimArea/dimLength/dimLength, Zero)
);

volScalarField chemistryhsSource
(
    IOobject
    (
        "chemistryhsSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("chemistryhsSource", dimEnergy/dimTime/dimVolume, 0.0)
);

dimensionedScalar pUnit("pUnit",dimensionSet (1,-1,-2,0,0,0,0),1.0);
dimensionedScalar TUnit("TUnit", dimTemperature, 1.0);
dimensionedScalar DUnit("DUnit",dimensionSet (0,2,-1,0,0,0,0),1.0);
dimensionedScalar lambdaUnit("lambdaUnit",dimensionSet (1,1,-3,-1,0,0,0),1.0);
dimensionedScalar MuUnit("MuUnit",dimensionSet (1,-1,-1,0,0,0,0),1.0);

PtrList<volScalarField> Diff(Y.size());
forAll(Y,i)
{
    word name("Diff"+Y[i].name());
    Diff.set
    (
        i,new volScalarField
        (
            IOobject
            (
                name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(name,dimensionSet (0,2,-1,0,0,0,0),0)
        )
    );
}

PtrList<volScalarField> d(0.5*Y.size()*(Y.size()+1));
forAll(d,i)
{
    word name("d"+i);
    d.set
    (
        i,new volScalarField
        (
            IOobject
            (
                name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(name,dimensionSet (0,0,0,0,0,0,0),0)
        )
    );
}

PtrList<volScalarField> X(Y.size());
forAll(X,i)
{
    word name("X"+Y[i].name());
    X.set
    (
        i,new volScalarField
        (
            IOobject
            (
                name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(name,dimensionSet (0,0,0,0,0,0,0),0)
        )
    );
}

//**********************Binary Diffusion****************//
List<dimensionedScalar> Diff1(0.5*Y.size()*(Y.size()+1));
List<dimensionedScalar> Diff2(0.5*Y.size()*(Y.size()+1));
List<dimensionedScalar> Diff3(0.5*Y.size()*(Y.size()+1));
List<dimensionedScalar> Diff4(0.5*Y.size()*(Y.size()+1));
//******************************************************//

volScalarField DiNumerator
(
    IOobject
    (
        "DiNumerator",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DiNumerator",dimensionSet (0,0,0,0,0,0,0),0.0)
);

volScalarField DiDenominator
(
    IOobject
    (
        "DiDenominator",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DiDenominator",dimensionSet (0,0,0,0,0,0,0),0.0)
);

volScalarField Di
(
    IOobject
    (
        "Di",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Di",dimensionSet (0,2,-1,0,0,0,0),0.0)
);

volScalarField Hsi
(
    IOobject
    (
        "Hsi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Hsi",dimensionSet (0,2,-2,0,0,0,0),0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField HEi
(
    IOobject
    (
        "HEi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Eei",dimensionSet (0,2,-2,0,0,0,0),0.0),
    zeroGradientFvPatchScalarField::typeName
);

volVectorField DiffError
(
    IOobject
    (
        "DiffError",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("DiffError",dimensionSet (1,-2,-1,0,0,0,0),vector::zero)
);

surfaceScalarField vk //volVectorField vk
(
    IOobject
    (
        "vk",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("vk",dimensionSet(1,2,-3,0,0,0,0),0.0)
);

volScalarField diffsh
(
    IOobject
    (
        "diffsh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("diffsh", dimEnergy/dimTime/dimVolume, 0.0)
);

//**********************thermal Diffusion****************//
List<dimensionedScalar> TDiff1(Y.size()*2);
List<dimensionedScalar> TDiff2(Y.size()*2);
List<dimensionedScalar> TDiff3(Y.size()*2);
List<dimensionedScalar> TDiff4(Y.size()*2);
//******************************************************//
volScalarField TDRatio_H
(
    IOobject
    (
        "TDRatio_H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TDRatio_H",dimensionSet (0,0,0,0,0,0,0),0.0)
);
volScalarField TDRatio_H2
(
    IOobject
    (
        "TDRatio_H2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("TDRatio_H2",dimensionSet (0,0,0,0,0,0,0),0.0)
);

//**********************lambda***********************//
PtrList<volScalarField> Lambda(Y.size());
forAll(Y,i)
{
    word name("Lambda"+Y[i].name());
    Lambda.set
    (
        i,new volScalarField
        (
            IOobject
            (
                name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(name,dimensionSet (1,1,-3,-1,0,0,0),0)
        )
    );
}

List<dimensionedScalar> Lambda1(Y.size());
List<dimensionedScalar> Lambda2(Y.size());
List<dimensionedScalar> Lambda3(Y.size());
List<dimensionedScalar> Lambda4(Y.size());

volScalarField mixLambdaA
(
    IOobject
    (
        "mixLambdaA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("mixLambdaA",dimensionSet (1,1,-3,-1,0,0,0),0.0)
);

volScalarField Lambdai
(
    IOobject
    (
        "Lambdai",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Lambdai",dimensionSet (1,1,-3,-1,0,0,0),0.0)
);

volScalarField mixLambdaB
(
    IOobject
    (
        "mixLambdaB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("mixLambdaB",dimensionSet (-1,-1,3,1,0,0,0),0.0)
);

volScalarField mixAlpha
(
    IOobject
    (
        "mixAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mixAlpha",dimensionSet (1,-1,-1,0,0,0,0),0.0)
);

//***********************viscosity*************************//
PtrList<volScalarField> Mu;
forAll(Y,i)
{
    word name("Mu"+Y[i].name());
    Mu.append
    (
        new volScalarField
        (
            IOobject
            (
                name,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(name,dimensionSet (1,-1,-1,0,0,0,0),0)
        )
    );
}

List<dimensionedScalar> Mu1(Y.size());
List<dimensionedScalar> Mu2(Y.size());
List<dimensionedScalar> Mu3(Y.size());
List<dimensionedScalar> Mu4(Y.size());

volScalarField mixMu
(
    IOobject
    (
        "mixMu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("mixMu",dimensionSet (1,-1,-1,0,0,0,0),0.0)
);

volScalarField yhi
(
    IOobject
    (
        "yhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("yhi",dimensionSet (0,0,0,0,0,0,0),0.0)
);

List<scalarList> groupSpecies(Y.size());
forAll(groupSpecies,i)
{
    groupSpecies[i].setSize(1,i);
}


if (mixtureAverageTransport)
{
    //***************** S. Zhong 20201008******************//
    IFstream speciesInfo(mesh.time().constant()/"trandat");
    if(!speciesInfo.good())
    {
        Info << FatalErrorIn( "cannot open trandat in constant") << abort(FatalError);
    }
    Foam::string   headerStr;
    Foam::word    thisSpecies;
    boolList foundSpecies(Y.size(), false);
    List<scalarList> trandata(Y.size());
    forAll(trandata,i)
    {
        trandata[i].setSize(6, 0.0);
    }

    speciesInfo.Foam::ISstream::getLine(headerStr);
    while( !speciesInfo.eof() )
    {
        speciesInfo >> thisSpecies;
        if(!speciesInfo.bad() )
        {
            for(label ii=0; ii< Y.size() ; ii++)
            {
                if(thisSpecies == Y[ii].name())
                {
                    foundSpecies[ii] = true;
                    speciesInfo >> trandata[ii][0] >> trandata[ii][1]  >> trandata[ii][2] >> trandata[ii][3]  >> trandata[ii][4]  >> trandata[ii][5]  ;
                    Info <<"read the trandata of "<< thisSpecies  << tab<<  trandata[ii][0] << tab << trandata[ii][1]<< tab << trandata[ii][2]<< tab << trandata[ii][3]<< tab<<trandata[ii][4]<<  tab<<trandata[ii][5]<< tab<< endl;
                }
            }
        }
        speciesInfo.Foam::ISstream::getLine(headerStr);
    }

    forAll(Y,i)
    {
        if(!foundSpecies[i])
        {
            Info << " species trandata for " << Y[i].name() << " is not found in constant/trandat "<<	FatalErrorIn( "something wrong in your constant/trandat") << abort(FatalError);
        }
    }
    // group the species that have the same trandata
    forAll(trandata,i)
    {
        for(label j = 0; j < i;j++)
        {
            if(i != j)
            {
                if(trandata[i] == trandata[j])
                {
                    label currentSize =  groupSpecies[j].size();
                    groupSpecies[j].resize(currentSize+1);
                    groupSpecies[j][currentSize] = i;
                    groupSpecies[i][0] = -1;
                    break;
                }
            }
        }
    }

    forAll(groupSpecies,i)
    {
        if(groupSpecies[i].size() > 1)
        {
            Info << "trandata group species:  " << endl;
            Info << "(" << endl;
            forAll(groupSpecies[i],j)
            {
                Info << Y[groupSpecies[i][j]].name() << endl;
            }
            Info << ")" << endl;
        }
    }

    Info << "trandata reduced species list:  " << endl;
    Info << "(" << endl;
    forAll(groupSpecies,i)
    {
        if(groupSpecies[i][0] != -1)
        {
            Info << Y[groupSpecies[i][0]].name() << endl;
        }
    }
    Info << ")" << endl;

    //*****************S. Zhong 2016.12.26******************//
    Info << "Read binaryDiff" << endl;
    IOdictionary binaryDiff
    (
        IOobject
        (
            "binaryDiff",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const dimensionedScalar dUnit
    (
        "dUnit",
        dimensionSet(0,0,0,0,0,0,0),
        1.0
    );

    boolList flags(0.5*Y.size()*(Y.size()+1), true);
    forAll(groupSpecies,i)
    {
        dictionary dDict;
        if(groupSpecies[i][0] != -1)
        {
            forAll(groupSpecies[i],h)
            {
                for(label j = i; j < groupSpecies.size();j ++)
                {
                    if(groupSpecies[j][0] != -1)
                    {
                        forAll(groupSpecies[j],z)
                        {
                            bool foundDict = false;
                            if(binaryDiff.found(Y[groupSpecies[i][h]].name()+"-"+Y[groupSpecies[j][z]].name()))
                            {
                                dDict = binaryDiff.subDict(Y[groupSpecies[i][h]].name()+"-"+Y[groupSpecies[j][z]].name());
                                foundDict = true;
                            }
                            else if(binaryDiff.found(Y[groupSpecies[j][z]].name()+"-"+Y[groupSpecies[i][h]].name()))
                            {
                                dDict = binaryDiff.subDict(Y[groupSpecies[j][z]].name()+"-"+Y[groupSpecies[i][h]].name());
                                foundDict = true;
                            }
                            else
                            {
                                label xx = z*h;
                                if(xx == (groupSpecies[j].size()*groupSpecies[i].size()))
                                {
                                    Info << "binary diffsusion coefficient for " << endl;
                                    Info << "(" << endl;
                                    forAll(groupSpecies[i],h)
                                    {
                                        Info << Y[groupSpecies[i][h]].name() << endl;
                                    }
                                    Info << ")" << endl;
                                    Info << "(" << endl;
                                    forAll(groupSpecies[j],z)
                                    {
                                        Info << Y[groupSpecies[j][z]].name() << endl;
                                    }
                                    Info << ")" << endl;
                                    Info << " are not provided in constant/binaryDiff" <<	FatalErrorIn( "something wrong in your constant/binaryDiff") << abort(FatalError);
                                }
                            }

                            if(foundDict)
                            {
                                forAll(groupSpecies[i],h)
                                {
                                    forAll(groupSpecies[j],z)
                                    {
                                        label m = max(groupSpecies[i][h],groupSpecies[j][z]);
                                        label n = min(groupSpecies[i][h],groupSpecies[j][z]);

                                        label k = Y.size()*n+m-0.5*n*(n+1);
                                        if(flags[k])
                                        {
                                            Diff1[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                                            Diff2[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                                            Diff3[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                                            Diff4[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

                                            Diff1[k] = readScalar(dDict.lookup("Diff1")) * dUnit;
                                            Diff2[k] = readScalar(dDict.lookup("Diff2")) * dUnit;
                                            Diff3[k] = readScalar(dDict.lookup("Diff3")) * dUnit;
                                            Diff4[k] = readScalar(dDict.lookup("Diff4")) * dUnit;
                                            flags[k] = false;
                                        }
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    Info << "species binary diffusion log polynomial fitting coefficient" << endl;
    Info << "(" << endl;
    forAll(Y,i)
    {
        for (label j = i; j < Y.size(); j++)
        {
            label m = max(i,j);
            label n = min(i,j);
            label k = Y.size()*n+m-0.5*n*(n+1);
            Info << Y[i].name() << "-" << Y[j].name() << tab << Diff1[k].value() << tab << Diff2[k].value() << tab << Diff3[k].value() << tab << Diff4[k].value() << endl;
        }
    }
    Info << ")" << endl;

    //********************S. Zhong 2017.3.14**************//
    Info << "Read speciesLambda" << endl;

    IOdictionary speciesLambda
    (
        IOobject
        (
            "speciesLambda",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const dimensionedScalar speciesLambdaUnit
    (
        "speciesLambdaUnit",
        dimensionSet(0,0,0,0,0,0,0),
        1.0
    );

    forAll(groupSpecies,i)
    {
        dictionary lambdaDict;
        if(groupSpecies[i][0] != -1)
        {
            forAll(groupSpecies[i],j)
            {
                if(speciesLambda.found(Y[groupSpecies[i][j]].name()))
                {
                    lambdaDict = speciesLambda.subDict(Y[groupSpecies[i][j]].name());
                    forAll(groupSpecies[i],j)
                    {
                        Lambda1[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                        Lambda2[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                        Lambda3[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                        Lambda4[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

                        Lambda1[groupSpecies[i][j]] = readScalar(lambdaDict.lookup("Lambda1")) * speciesLambdaUnit;
                        Lambda2[groupSpecies[i][j]] = readScalar(lambdaDict.lookup("Lambda2")) * speciesLambdaUnit;
                        Lambda3[groupSpecies[i][j]] = readScalar(lambdaDict.lookup("Lambda3")) * speciesLambdaUnit;
                        Lambda4[groupSpecies[i][j]] = readScalar(lambdaDict.lookup("Lambda4")) * speciesLambdaUnit;
                    }
                    break;
                }
                else
                {
                    if( j == (groupSpecies[i].size()-1) )
                    {
                        Info << "species conductivity log polynomial fitting coefficient for " << endl;
                        Info << "(" << endl;
                        forAll(groupSpecies[i],j)
                        {
                            Info << Y[groupSpecies[i][j]].name() << endl;
                        }
                        Info << ")" << endl;
                        Info << " are not provided in constant/speciesLambda" <<	FatalErrorIn( "something wrong in your constant/speciesLambda") << abort(FatalError);
                    }
                }
            }
        }
    }

    Info << "species conductivity log polynomial fitting coefficient" << endl;
    Info << "(" << endl;
    forAll(Y,i)
    {
        Info << Y[i].name() << tab << Lambda1[i].value() << tab << Lambda2[i].value() << tab << Lambda3[i].value() << tab << Lambda4[i].value() << endl;
    }
    Info << ")" << endl;

    //********************S. Zhong 2017.3.14**************//
    Info << "Read speciesMu" << endl;

    IOdictionary speciesMu
    (
        IOobject
        (
            "speciesMu",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const dimensionedScalar speciesMuUnit
    (
        "speciesMuUnit",
        dimensionSet(0,0,0,0,0,0,0),
        1.0
    );

    forAll(groupSpecies,i)
    {
        dictionary muDict;
        if(groupSpecies[i][0] != -1)
        {
            forAll(groupSpecies[i],j)
            {
                if(speciesMu.found(Y[groupSpecies[i][j]].name()))
                {
                    muDict = speciesMu.subDict(Y[groupSpecies[i][j]].name());
                    forAll(groupSpecies[i],j)
                    {
                        Mu1[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                        Mu2[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                        Mu3[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                        Mu4[groupSpecies[i][j]].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

                        Mu1[groupSpecies[i][j]] = readScalar(muDict.lookup("Mu1")) * speciesMuUnit;
                        Mu2[groupSpecies[i][j]] = readScalar(muDict.lookup("Mu2")) * speciesMuUnit;
                        Mu3[groupSpecies[i][j]] = readScalar(muDict.lookup("Mu3")) * speciesMuUnit;
                        Mu4[groupSpecies[i][j]] = readScalar(muDict.lookup("Mu4")) * speciesMuUnit;
                    }
                    break;
                }
                else
                {
                    if( j == (groupSpecies[i].size()-1) )
                    {
                        Info << "species viscosity log polynomial fitting coefficient for " << endl;
                        Info << "(" << endl;
                        forAll(groupSpecies[i],j)
                        {
                            Info << Y[groupSpecies[i][j]].name() << endl;
                        }
                        Info << ")" << endl;
                        Info << " are not provided in constant/speciesMu" <<	FatalErrorIn( "something wrong in your constant/speciesMu") << abort(FatalError);
                    }
                }
            }
        }
    }

    Info << "species viscosity log polynomial fitting coefficient" << endl;
    Info << "(" << endl;
    forAll(Y,i)
    {
        Info << Y[i].name() << tab << Mu1[i].value() << tab << Mu2[i].value() << tab << Mu3[i].value() << tab << Mu4[i].value() << endl;
    }
    Info << ")" << endl;

    //*****************S. Zhong 2020.09.08******************//
    Info << "Read thermoDiff" << endl;
    IOdictionary thermoDiff
    (
        IOobject
        (
            "thermoDiff",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const dimensionedScalar tdUnit
    (
        "tdUnit",
        dimensionSet(0,0,0,0,0,0,0),
        1.0
    );
    boolList tflags(2*Y.size(), true);
    forAll(Y,i)
    {
        if(groupSpecies[i][0] != -1)
        {
            if(Y[i].name() == "H" || Y[i].name() == "H2")
            {
                dictionary tdDict;
                forAll(groupSpecies,j)
                {
                    if(groupSpecies[j][0] != -1)
                    {
                        forAll(groupSpecies[j],z)
                        {
                            bool foundDict = false;
                            if(thermoDiff.found(Y[i].name()+"-"+Y[groupSpecies[j][z]].name()))
                            {
                                tdDict = thermoDiff.subDict(Y[i].name()+"-"+Y[groupSpecies[j][z]].name());
                                foundDict = true;
                            }
                            else if(thermoDiff.found(Y[groupSpecies[j][z]].name()+"-"+Y[i].name()))
                            {
                                tdDict = thermoDiff.subDict(Y[groupSpecies[j][z]].name()+"-"+Y[i].name());
                                foundDict = true;
                            }
                            else
                            {
                                if( z == (groupSpecies[j].size()-1) )
                                {
                                    Info << "thermal diffsusion coefficient for " << Y[i].name() << " and " << endl;
                                    Info << "(" << endl;
                                    forAll(groupSpecies[j],z)
                                    {
                                        Info << Y[groupSpecies[j][z]].name() << endl;
                                    }
                                    Info << ")" << endl;
                                    Info << " are not provided in constant/thermoDiff" <<	FatalErrorIn( "something wrong in your constant/thermoDiff") << abort(FatalError);
                                }
                            }

                            if(foundDict)
                            {
                                forAll(groupSpecies[j],z)
                                {
                                    if( i != groupSpecies[j][z])
                                    {
                                        label k = -1;
                                        if(Y[i].name() == "H")
                                        {
                                            k = groupSpecies[j][z];
                                        }
                                        if(Y[i].name() == "H2")
                                        {
                                            k = Y.size()+groupSpecies[j][z];
                                        }
                                        if(tflags[k])
                                        {
                                            TDiff1[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                                            TDiff2[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                                            TDiff3[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));
                                            TDiff4[k].dimensions().reset(dimensionSet(0,0,0,0,0,0,0));

                                            TDiff1[k] = readScalar(tdDict.lookup("ThermDiff_1")) * tdUnit;
                                            TDiff2[k] = readScalar(tdDict.lookup("ThermDiff_2")) * tdUnit;
                                            TDiff3[k] = readScalar(tdDict.lookup("ThermDiff_3")) * tdUnit;
                                            TDiff4[k] = readScalar(tdDict.lookup("ThermDiff_4")) * tdUnit;
                                            tflags[k] = false;
                                        }
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    Info << "species (H and H2, molecular weight less than 5) thermal diffusion ratio polynomial fitting coefficient" << endl;
    Info << "(" << endl;
    forAll(Y,i)
    {
        if(Y[i].name() == "H" || Y[i].name() == "H2")
        {
            for (label j = 0; j < Y.size(); j++)
            {
                label k = -1;
                if(Y[i].name() == "H")
                {
                    k = j;
                    Info << Y[i].name() << "-" << Y[j].name() << tab << TDiff1[k].value() << tab << TDiff2[k].value() << tab << TDiff3[k].value() << tab << TDiff4[k].value() << endl;
                }
                if(Y[i].name() == "H2")
                {
                    k = Y.size()+j;
                    Info << Y[i].name() << "-" << Y[j].name() << tab << TDiff1[k].value() << tab << TDiff2[k].value() << tab << TDiff3[k].value() << tab << TDiff4[k].value() << endl;
                }
            }
        }
    }
    Info << ")" << endl;
}

