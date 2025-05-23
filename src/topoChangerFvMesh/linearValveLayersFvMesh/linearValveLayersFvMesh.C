/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "linearValveLayersFvMesh.H"
#include "Time.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "pointField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearValveLayersFvMesh, 0);
    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        linearValveLayersFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::linearValveLayersFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if
    (
        pointZones().size()
     || faceZones().size()
     || cellZones().size()
     || topoChanger_.size()
    )
    {
        InfoInFunction
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    // Add zones
    List<pointZone*> pz(1);
    List<faceZone*> fz(4);
    List<cellZone*> cz(0);


    // An empty zone for cut points
    pz[0] = new pointZone("cutPointZone", 0, pointZones());


    // Do face zones for slider

    // Inner slider
    const word innerSliderName
    (
        motionDict_.subDict("slider").get<word>("inside")
    );
    const polyPatch& innerSlider = boundaryMesh()[innerSliderName];

    fz[0] = new faceZone
    (
        "insideSliderZone",
        identity(innerSlider.range()),
        false, // none are flipped
        0,
        faceZones()
    );

    // Outer slider
    const word outerSliderName
    (
        motionDict_.subDict("slider").get<word>("outside")
    );
    const polyPatch& outerSlider = boundaryMesh()[outerSliderName];

    fz[1] = new faceZone
    (
        "outsideSliderZone",
        identity(outsideSlider.range()),
        false, // none are flipped
        1,
        faceZones()
    );

    // An empty zone for cut faces
    fz[2] = new faceZone("cutFaceZone", 2, faceZones());

    // Add face zone for layer addition
    const word layerPatchName
    (
        motionDict_.subDict("layer").get<word>("patch")
    );

    const polyPatch& layerPatch = boundaryMesh()[layerPatchName];

    fz[3] = new faceZone
    (
        "valveLayerZone",
        identity(layerPatch.range()),
        lpf,
        true, // all are flipped
        0,
        faceZones()
    );


    Info<< "Adding point and face zones" << endl;
    addZones(pz, fz, cz);

    // Add a topology modifier

    List<polyMeshModifier*> tm(2);

    tm[0] = new slidingInterface
    (
        "valveSlider",
        0,
        topoChanger_,
        outerSliderName + "Zone",
        innerSliderName + "Zone",
        "cutPointZone",
        "cutFaceZone",
        outerSliderName,
        innerSliderName,
        slidingInterface::INTEGRAL,
        true                          // Attach-detach action
    );

    tm[1] =
        new layerAdditionRemoval
        (
            "valveLayer",
            1,
            topoChanger_,
            "valveLayerZone",
            motionDict_.subDict("layer").get<scalar>("minThickness"),
            motionDict_.subDict("layer").get<scalar>("maxThickness")
        );


    Info<< "Adding topology modifiers" << endl;
    addTopologyModifiers(tm);

    // Write mesh
    write();
}


void Foam::linearValveLayersFvMesh::makeLayersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable layering
    forAll(topoChanges, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].disable();
        }
        else
        {
            FatalErrorInFunction
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::linearValveLayersFvMesh::makeSlidersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable sliding interface
    forAll(topoChanges, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanges[modI]))
        {
            topoChanges[modI].disable();
        }
        else if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else
        {
            FatalErrorInFunction
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


bool Foam::linearValveLayersFvMesh::attached() const
{
    const polyTopoChanger& topoChanges = topoChanger_;

    bool result = false;

    forAll(topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            result =
                result
             || refCast<const slidingInterface>(topoChanges[modI]).attached();
        }
    }

    // Check thal all sliders are in sync (debug only)
    forAll(topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            if
            (
                result
             != refCast<const slidingInterface>(topoChanges[modI]).attached()
            )
            {
                FatalErrorInFunction
                    << "Slider " << modI << " named "
                    << topoChanges[modI].name()
                    << " out of sync: Should be" << result
                    << abort(FatalError);
            }
        }
    }

    return result;
}


Foam::tmp<Foam::pointField> Foam::linearValveLayersFvMesh::newPoints() const
{
    auto tnewPoints = tmp<pointField>::New(points());
    auto& np = tnewPoints();

    const word layerPatchName
    (
        motionDict_.subDict("layer").get<word>("patch")
    );

    const polyPatch& layerPatch = boundaryMesh()[layerPatchName];

    const labelList& patchPoints = layerPatch.meshPoints();

    const vector vel
    (
        motionDict_.get<vector>("pistonVelocity")
    );

    forAll(patchPoints, ppI)
    {
        np[patchPoints[ppI]] += vel*time().deltaTValue();
    }

    return tnewPoints;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::linearValveLayersFvMesh::linearValveLayersFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io),
    motionDict_
    (
        IOdictionary::readContents
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ
            )
        ).optionalSubDict(typeName + "Coeffs")
    )
{
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearValveLayersFvMesh::~linearValveLayersFvMesh()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::linearValveLayersFvMesh::update()
{
    // Detaching the interface
    if (attached())
    {
        Info<< "Decoupling sliding interfaces" << endl;
        makeSlidersLive();

        // Changing topology
        resetMorph();
        setMorphTimeIndex(3*time().timeIndex());
        updateMesh();
    }
    else
    {
        Info<< "Sliding interfaces decoupled" << endl;
    }

    // Perform layer action and mesh motion
    makeLayersLive();

    // Changing topology
    resetMorph();
    setMorphTimeIndex(3*time().timeIndex() + 1);
    updateMesh();

    if (topoChangeMap)
    {
        if (topoChangeMap().hasMotionPoints())
        {
            Info<< "Topology change; executing pre-motion" << endl;
            movePoints(topoChangeMap().preMotionPoints());
        }
    }

    // Move points
    movePoints(newPoints());

    // Attach the interface
    Info<< "Coupling sliding interfaces" << endl;
    makeSlidersLive();

    // Changing topology
    resetMorph();
    setMorphTimeIndex(3*time().timeIndex() + 2);
    updateMesh();

    //Info<< "Moving points post slider attach" << endl;
    //const pointField p = allPoints();
    //movePoints(p);

    Info<< "Sliding interfaces coupled: " << attached() << endl;
}


// ************************************************************************* //
