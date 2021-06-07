/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "forces.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forces, 0);
    addToRunTimeSelectionTable(functionObject, forces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::word Foam::functionObjects::forces::fieldName(const word& name) const
{
    return this->name() + ":" + name;
}


void Foam::functionObjects::forces::setCoordinateSystem
(
    const dictionary& dict,
    const word& e3Name,
    const word& e1Name
)
{
    coordSys_.clear();

    if (dict.readIfPresent<point>("CofR", coordSys_.origin()))
    {
        const vector e3 = e3Name == word::null ?
            vector(0, 0, 1) : dict.get<vector>(e3Name);
        const vector e1 = e1Name == word::null ?
            vector(1, 0, 0) : dict.get<vector>(e1Name);

        coordSys_ = coordSystem::cartesian(coordSys_.origin(), e3, e1);
    }
    else
    {
        // The 'coordinateSystem' sub-dictionary is optional,
        // but enforce use of a cartesian system if not found.

        if (dict.found(coordinateSystem::typeName_()))
        {
            // New() for access to indirect (global) coordinate system
            coordSys_ =
                coordinateSystem::New
                (
                    obr_,
                    dict,
                    coordinateSystem::typeName_()
                );
        }
        else
        {
            coordSys_ = coordSystem::cartesian(dict);
        }
    }
}


void Foam::functionObjects::forces::initialise()
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database"
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            !foundObject<volVectorField>(UName_)
         || !foundObject<volScalarField>(pName_)
        )
        {
            FatalErrorInFunction
                << "Could not find U: " << UName_
                << " or p:" << pName_ << " in database"
                << exit(FatalError);
        }

        if (rhoName_ != "rhoInf" && !foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Could not find rho:" << rhoName_ << " in database"
                << exit(FatalError);
        }
    }

    binModelPtr_->initialise();
    const label nBin = binModelPtr_->nBin();

    // Set size of force and moment containers according to the bin
    for (direction i = 0; i < vector::nComponents; ++i)
    {
        forces_[i].setSize(nBin, Zero);
        moments_[i].setSize(nBin, Zero);
    }

    initialised_ = true;
}


void Foam::functionObjects::forces::reset()
{
    forces_ = Zero;
    moments_ = Zero;

    if (writeFields_)
    {
        auto& force = lookupObjectRef<volVectorField>(fieldName("force"));

        force == dimensionedVector(force.dimensions(), Zero);

        auto& moment = lookupObjectRef<volVectorField>(fieldName("moment"));

        moment == dimensionedVector(moment.dimensions(), Zero);
    }
}


Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::forces::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    const auto& U = lookupObject<volVectorField>(UName_);

    if (foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff(U);
    }
    else if (foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return rho()*turb.devReff(U);
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo = lookupObject<fluidThermo>(fluidThermo::dictName);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
            lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        return -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::mu() const
{
    if (foundObject<fluidThermo>(basicThermo::dictName))
    {
        const auto& thermo = lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
             lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::forces::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                "rho",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimDensity, rhoRef_)
        );
    }

    return (lookupObject<volScalarField>(rhoName_));
}


Foam::scalar Foam::functionObjects::forces::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1;
    }

    if (rhoName_ != "rhoInf")
    {
        FatalErrorInFunction
            << "Dynamic pressure is expected but kinematic is provided."
            << exit(FatalError);
    }

    return rhoRef_;
}


void Foam::functionObjects::forces::addToFields
(
    const label patchi,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT
)
{
    if (!writeFields_)
    {
        return;
    }

    auto& force = lookupObjectRef<volVectorField>(fieldName("force"));
    vectorField& pf = force.boundaryFieldRef()[patchi];
    pf += fN + fT;

    auto& moment = lookupObjectRef<volVectorField>(fieldName("moment"));
    vectorField& pm = moment.boundaryFieldRef()[patchi];
    pm = Md^pf;
}


void Foam::functionObjects::forces::addToFields
(
    const labelList& cellIDs,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    if (!writeFields_)
    {
        return;
    }

    auto& force = lookupObjectRef<volVectorField>(fieldName("force"));
    auto& moment = lookupObjectRef<volVectorField>(fieldName("moment"));

    forAll(cellIDs, i)
    {
        const label celli = cellIDs[i];
        force[celli] += fN[i] + fT[i] + fP[i];
        moment[celli] = Md[i]^force[celli];
    }
}


Foam::List<Foam::vector> Foam::functionObjects::forces::integrateData
(
    const vectorField& fm0,
    const vectorField& fm1,
    const vectorField& fm2
) const
{
    vector pressure(sum(fm0));
    vector viscous(sum(fm1));
    vector porous = porosity_ ? sum(fm2) : Zero;
    vector total(pressure + viscous + porous);

    return List<vector>({total, pressure, viscous, porous});
}


bool Foam::functionObjects::forces::createDataFiles()
{
    if (!forceFilePtr_.valid())
    {
        forceFilePtr_ = createFile("force");
        writeDataFileHeader("Force", forceFilePtr_());
    }

    if (!momentFilePtr_.valid())
    {
        momentFilePtr_ = createFile("moment");
        writeDataFileHeader("Moment", momentFilePtr_());
    }

    return true;
}


bool Foam::functionObjects::forces::writeDataFileHeader
(
    const word& header,
    OFstream& os
) const
{
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", coordSys_.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "(total_x total_y total_z)");
    writeTabbed(os, "(pressure_x pressure_y pressure_z)");
    writeTabbed(os, "(viscous_x viscous_y viscous_z)");

    if (porosity_)
    {
        writeTabbed(os, "(porous_x porous_y porous_z)");
    }

    os  << endl;

    return true;
}


bool Foam::functionObjects::forces::writeDataFiles()
{
    writeDataFile(sumForces_, forceFilePtr_);

    writeDataFile(sumMoments_, momentFilePtr_);

    return true;
}


bool Foam::functionObjects::forces::writeDataFile
(
    const List<vector>& fm,
    autoPtr<OFstream>& osPtr
) const
{
    OFstream& os = osPtr();

    writeCurrentTime(os);

    os  << tab << fm[0]
        << tab << fm[1]
        << tab << fm[2];

    if (porosity_)
    {
        os  << tab << fm[3];
    }

    os  << endl;

    return true;
}


void Foam::functionObjects::forces::logData
(
    const string& descriptor,
    const List<vector>& fm
) const
{
    if (!log)
    {
        return;
    }

    Log << "    Sum of " << descriptor.c_str() << nl
        << "        Total    : " << fm[0] << nl
        << "        Pressure : " << fm[1] << nl
        << "        Viscous  : " << fm[2] << nl;

    if (porosity_)
    {
        Log << "        Porous   : " << fm[3] << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forces::forces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    binModelPtr_(binModel::New(name, dict, mesh_)),
    forces_(vector::nComponents),
    moments_(vector::nComponents),
    sumForces_(Zero),
    sumMoments_(Zero),
    forceFilePtr_(),
    momentFilePtr_(),
    coordSys_(),
    patchSet_(),
    rhoRef_(VGREAT),
    pRef_(0),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    fDName_("fD"),
    directForceDensity_(false),
    porosity_(false),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


Foam::functionObjects::forces::forces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),
    binModelPtr_(binModel::New(name, dict, mesh_)),
    forces_(vector::nComponents),
    moments_(vector::nComponents),
    sumForces_(Zero),
    sumMoments_(Zero),
    forceFilePtr_(),
    momentFilePtr_(),
    coordSys_(),
    patchSet_(),
    rhoRef_(VGREAT),
    pRef_(0),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    fDName_("fD"),
    directForceDensity_(false),
    porosity_(false),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forces::read(const dictionary& dict)
{
    if (!(fvMeshFunctionObject::read(dict) && writeFile::read(dict)))
    {
        return false;
    }

    initialised_ = false;

    Info<< type() << " " << name() << ":" << endl;

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        );

    dict.readIfPresent("directForceDensity", directForceDensity_);
    if (directForceDensity_)
    {
        // Optional name entry for fD
        if (dict.readIfPresent<word>("fD", fDName_))
        {
            Info<< "    fD: " << fDName_ << endl;
        }
    }
    else
    {
        // Optional field name entries
        if (dict.readIfPresent<word>("p", pName_))
        {
            Info<< "    p: " << pName_ << endl;
        }
        if (dict.readIfPresent<word>("U", UName_))
        {
            Info<< "    U: " << UName_ << endl;
        }
        if (dict.readIfPresent<word>("rho", rhoName_))
        {
            Info<< "    rho: " << rhoName_ << endl;
        }

        // Reference density needed for incompressible calculations
        if (rhoName_ == "rhoInf")
        {
            rhoRef_ = dict.getCheck<scalar>("rhoInf", scalarMinMax::ge(SMALL));
            Info<< "    Freestream density (rhoInf) set to " << rhoRef_ << endl;
        }

        // Reference pressure, 0 by default
        if (dict.readIfPresent<scalar>("pRef", pRef_))
        {
            Info<< "    Reference pressure (pRef) set to " << pRef_ << endl;
        }
    }

    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        Info<< "    Including porosity effects" << endl;
    }
    else
    {
        Info<< "    Not including porosity effects" << endl;
    }

    writeFields_ = dict.getOrDefault("writeFields", false);
    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;

        auto* forcePtr
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("force"),
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector(dimForce, Zero)
            )
        );

        mesh_.objectRegistry::store(forcePtr);

        auto* momentPtr
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("moment"),
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector(dimForce*dimLength, Zero)
            )
        );

        mesh_.objectRegistry::store(momentPtr);
    }

    binModelPtr_->read(dict);

    return true;
}


void Foam::functionObjects::forces::calcForcesMoments()
{
    initialise();

    reset();

    if (directForceDensity_)
    {
        const auto& fD = lookupObject<volVectorField>(fDName_);

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        for (const label patchi : patchSet_)
        {
            const vectorField& d = mesh_.C().boundaryField()[patchi];

            const vectorField Md(d - coordSys_.origin());

            const scalarField sA(mag(Sfb[patchi]));

            // Normal force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            const vectorField fN
            (
                Sfb[patchi]/sA
               *(
                    Sfb[patchi] & fD.boundaryField()[patchi]
                )
            );

            // Tangential force (total force minus normal fN)
            const vectorField fT(sA*fD.boundaryField()[patchi] - fN);

            addToFields(patchi, Md, fN, fT);

            binModelPtr_->applyBins(forces_, moments_, d, Md, fN, fT);
        }
    }
    else
    {
        const auto& p = lookupObject<volScalarField>(pName_);

        const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

        tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
        const volSymmTensorField::Boundary& devRhoReffb
            = tdevRhoReff().boundaryField();

        // Scale pRef by density for incompressible simulations
        const scalar rhoRef = rho(p);
        const scalar pRef = pRef_/rhoRef;

        for (const label patchi : patchSet_)
        {
            const vectorField& d = mesh_.C().boundaryField()[patchi];

            const vectorField Md(d - coordSys_.origin());

            const vectorField fN
            (
                rhoRef*Sfb[patchi]*(p.boundaryField()[patchi] - pRef)
            );

            const vectorField fT(Sfb[patchi] & devRhoReffb[patchi]);

            addToFields(patchi, Md, fN, fT);

            binModelPtr_->applyBins(forces_, moments_, d, Md, fN, fT);
        }
    }

    if (porosity_)
    {
        const auto& U = lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, "
                << "but no porosity models found in the database"
                << endl;
        }

        forAllConstIters(models, iter)
        {
            // Non-const access required if mesh is changing
            auto& pm = const_cast<porosityModel&>(*iter());

            const vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            for (const label zonei : cellZoneIDs)
            {
                const cellZone& cZone = mesh_.cellZones()[zonei];

                const vectorField d(mesh_.C(), cZone);
                const vectorField fP(fPTot, cZone);
                const vectorField Md(d - coordSys_.origin());

                const vectorField fNT(Md.size(), Zero);

                addToFields(cZone, Md, fNT, fNT, fP);

                binModelPtr_->applyBins(forces_, moments_, d, Md, fNT, fNT, fP);
            }
        }
    }

    Pstream::listCombineGather(forces_, plusEqOp<vectorField>());
    Pstream::listCombineGather(moments_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(forces_);
    Pstream::listCombineScatter(moments_);
}


Foam::vector Foam::functionObjects::forces::forceEff() const
{
    return sum(forces_[0]) + sum(forces_[1]) + sum(forces_[2]);
}


Foam::vector Foam::functionObjects::forces::momentEff() const
{
    return sum(moments_[0]) + sum(moments_[1]) + sum(moments_[2]);
}


bool Foam::functionObjects::forces::execute()
{
    calcForcesMoments();

    Log << type() << " " << name() << " write:" << nl;

    sumForces_ =
        integrateData
        (
            coordSys_.localVector(forces_[0]),
            coordSys_.localVector(forces_[1]),
            coordSys_.localVector(forces_[2])
        );

    logData("forces", sumForces_);

    sumMoments_ =
        integrateData
        (
            coordSys_.localVector(moments_[0]),
            coordSys_.localVector(moments_[1]),
            coordSys_.localVector(moments_[2])
        );

    logData("moments", sumMoments_);

    // Register forces and moments data in the local coordinate system
    const List<word> dirs({"normal", "tangential", "porous"});
    for (label i = 0; i < vector::nComponents; ++i)
    {
        setResult(dirs[i] + "Force", sumForces_[i+1]);
        setResult(dirs[i] + "Moment", sumMoments_[i+1]);
    }

    return true;
}


bool Foam::functionObjects::forces::write()
{
    if (writeToFile() && Pstream::master())
    {
        Log << "    writing force and moment files." << endl;

        createDataFiles();

        binModelPtr_->createBinnedDataFiles();

        writeDataFiles();

        binModelPtr_->writeBinnedDataFiles(forces_, moments_, coordSys_);

        Log << endl;
    }

    return true;
}


// ************************************************************************* //
