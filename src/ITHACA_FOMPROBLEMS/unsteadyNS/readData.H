Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

scalar density =
dimensionedScalar
(
    "density",
    dimless,
    transportProperties.lookup("density")
).value();


scalar Rtot =
dimensionedScalar
(
    "Rtot",
    dimless,
    transportProperties.lookup("Rtot")
).value();

scalar Tau =
dimensionedScalar
(
    "Tau",
    dimless,
    transportProperties.lookup("Tau")
).value();


scalar Ctot =
dimensionedScalar
(
    "Ctot",
    dimless,
    transportProperties.lookup("Ctot")
).value();



/*Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Info<< "Reading field Uav\n" << endl;
volVectorField Uav
(
    IOobject
    (
        "Uav",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("Uav", dimensionSet(0,1,-1,0,0,0,0), vector::zero)

);

//Info<< "Reading field WSS\n" << endl;
volVectorField WSS
(
    IOobject
    (
        "WSS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("WSS", dimensionSet(0,2,-2,0,0,0,0), vector::zero)
);

//Info<< "Reading field TAWSS\n" << endl;
volScalarField TAWSS
(
    IOobject
    (
        "TAWSS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("TAWSS", dimensionSet(0,2,-2,0,0,0,0), 0.0)
);*/

/*volVectorField TAWSS_2
(
    IOobject
    (
        "TAWSS_2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("TAWSS_2", dimensionSet(0,2,-2,0,0,0,0), vector::zero)
);

volScalarField OSI
(
    IOobject
    (
        "OSI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("OSI", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

volScalarField RRT
(
    IOobject
    (
        "RRT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("RRT", dimensionSet(0,-2,2,0,0,0,0), 0.0)
);


volScalarField LNH
(
    IOobject
    (
        "LNH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("LNH", dimensionSet(0,0,0,0,0,0,0), 0.0)
);


volScalarField LNH_av
(
    IOobject
    (
        "LNH_av",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("LNH_av", dimensionSet(0,0,0,0,0,0,0), 0.0)
);
*/

//#include "createPhi.H"


/*label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());*/
