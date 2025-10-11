/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#include "regionFaModel.H"
#include "faMesh.H"
#include "faMeshesRegistry.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(regionFaModel, 0);
}
}

const Foam::word
Foam::regionModels::regionFaModel::regionFaModelName("regionFaModel");


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Return IOobject with name qualified with region and area names
Foam::IOobject createModelIOobject
(
    const Foam::polyMesh& mesh,
    // const Foam::word& baseName, <- always regionFaModelName
    const Foam::word& regionName,
    const Foam::word& areaName
)
{
    using namespace Foam;

    // Default: regionFaModel.<regionName>
    word objName = IOobject::groupName
    (
        Foam::regionModels::regionFaModel::regionFaModelName,
        regionName
    );

    // Append '.<area-name>' or nothing
    objName.ext(polyMesh::regionName(areaName));

    return IOobject
    (
        objName,
        mesh.time().constant(),
        faMeshesRegistry::New(mesh).thisDb(),
        IOobjectOption::NO_READ,
        IOobjectOption::NO_WRITE,
        IOobjectOption::REGISTER
    );
}


// Return IOobject with name qualified with region and area names
Foam::IOobject createPropertiesIOobject
(
    const Foam::polyMesh& mesh,
    // const Foam::word& baseName, <- always regionFaModelName
    const Foam::word& regionName,
    const Foam::word& areaName
)
{
    using namespace Foam;

    const fileName uniformPath
    (
        word("uniform")
      / Foam::regionModels::regionFaModel::regionFaModelName
    );

    const word objName
    (
        IOobject::groupName
        (
            (regionName + "OutputProperties"),
            polyMesh::regionName(areaName)
        )
    );

    // NOTE (2025-10-01):
    // Cannot hold the OutputProperties within
    //    - faMeshesRegistry::New(mesh).thisDb()
    // since this produces a uniform path that we do not yet handle
    //
    // ->    "<time>/finite-area/uniform/regionFaModel/<model-region>"
    // vs:   "<time>/uniform/regionFaModel/<model-region>"
    //
    // The difference being that we only look for 'uniform' at the
    // first sub-level within the time directory when decomposing etc.

    IOobject legacy
    (
        objName,
        mesh.time().timeName(),
        uniformPath/regionName,

        // Not possible:  faMeshesRegistry::New(mesh).thisDb(),
        mesh,  // Registered on volume mesh!

        IOobjectOption::READ_IF_PRESENT,
        IOobjectOption::NO_WRITE,
        IOobjectOption::REGISTER
    );

    return legacy;
}

} // End anonymous namespace


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::regionModels::regionFaModel::constructMeshObjects()
{
    regionMeshPtr_.reset
    (
        new faMesh(areaName_, primaryMesh_)
    );
}


void Foam::regionModels::regionFaModel::initialise()
{
    if (debug)
    {
        Pout<< "regionFaModel::initialise()" << endl;
    }

    vsmPtr_.reset(new volSurfaceMapping(regionMeshPtr_()));

    if (!outputPropertiesPtr_)
    {
        outputPropertiesPtr_.reset
        (
            new IOdictionary
            (
                createPropertiesIOobject
                (
                    primaryMesh_,
                    // regionFaModelName,
                    regionName_,
                    areaName_
                )
            )
        );
    }
}


bool Foam::regionModels::regionFaModel::init(const dictionary& dict)
{
    if (active_)
    {
        if (const dictionary* dictptr = dict.findDict(modelName_ + "Coeffs"))
        {
            coeffs_ <<= *dictptr;
        }

        infoOutput_.readIfPresent("infoOutput", dict);

        return true;
    }

    return false;
}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

const Foam::volSurfaceMapping& Foam::regionModels::regionFaModel::vsm() const
{
    return vsmPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::regionFaModel::regionFaModel
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    const dictionary& dict,
    bool readFields
)
:
    IOdictionary
    (
        createModelIOobject
        (
            mesh,
            // regionFaModelName,
            dict.get<word>("region"),
            dict.getOrDefault<word>("area", polyMesh::defaultRegion)
        )
    ),
    primaryMesh_(mesh),
    time_(mesh.time()),
    active_(dict.get<Switch>("active")),
    infoOutput_(false),
    modelName_(modelName),
    areaName_(dict.getOrDefault<word>("area", polyMesh::defaultRegion)),
    regionName_(dict.get<word>("region")),
    coeffs_(dict.subOrEmptyDict(modelName + "Coeffs"))
{
    // Suffix hint for variable names
    if
    (
        coeffs_.readIfPresent("suffixing", suffixHint_)
     || dict.readIfPresent("suffixing", suffixHint_)
    )
    {
        Switch sw = Switch::find(suffixHint_);

        if (sw.good())
        {
            if (!sw)  // No suffix
            {
                suffixHint_.clear();
            }
        }
        else if (suffixHint_ == "default")
        {
            // Use default suffix
            sw = true;
        }

        if (sw)  // Default suffix
        {
            suffixHint_ = '_' + regionName_;
        }
    }

    constructMeshObjects();
    initialise();

    if (readFields)
    {
        init(dict);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::regionModels::regionFaModel::evolve()
{
    if (active_)
    {
        Info<< "\nEvolving " << modelName_
            << " for region " << regionMesh().name();

        if (!polyMesh::regionName(areaName_).empty())
        {
            Info<< " [" << areaName_ << "]";
        }
        Info<< endl;

        preEvolveRegion();

        evolveRegion();

        postEvolveRegion();

        // Provide some feedback
        if (infoOutput_)
        {
            Info<< incrIndent;
            info();
            Info<< decrIndent << endl;
        }
    }
}


void Foam::regionModels::regionFaModel::preEvolveRegion()
{}


void Foam::regionModels::regionFaModel::evolveRegion()
{}


void Foam::regionModels::regionFaModel::postEvolveRegion()
{}


Foam::scalar Foam::regionModels::regionFaModel::CourantNumber() const
{
    return 0;
}


// ************************************************************************* //
