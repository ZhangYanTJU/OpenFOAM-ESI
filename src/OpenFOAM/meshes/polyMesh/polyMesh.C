/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017, 2020 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "polyMesh.H"
#include "Time.H"
#include "cellIOList.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "globalMeshData.H"
#include "processorPolyPatch.H"
#include "polyMeshTetDecomposition.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "MeshObject.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyMesh, 0);
}

Foam::word Foam::polyMesh::defaultRegion = "region0";

Foam::word Foam::polyMesh::meshSubDir = "polyMesh";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMesh::calcDirections() const
{
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        solutionD_[cmpt] = 1;
    }

    // Knock out empty and wedge directions. Note:they will be present on all
    // domains.

    bool hasEmptyPatches = false;
    bool hasWedgePatches = false;

    vector emptyDirVec = Zero;
    vector wedgeDirVec = Zero;

    forAll(boundaryMesh(), patchi)
    {
        const polyPatch& pp = boundaryMesh()[patchi];
        if (isA<emptyPolyPatch>(pp))
        {
            // Force calculation of geometric properties, independent of
            // size. This avoids parallel synchronisation problems.
            const vectorField::subField fa(pp.faceAreas());

            if (pp.size())
            {
                hasEmptyPatches = true;
                emptyDirVec += sum(cmptMag(fa));
            }
        }
        else if (isA<wedgePolyPatch>(pp))
        {
            const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>(pp);

            // Force calculation of geometric properties, independent of
            // size. This avoids parallel synchronisation problems.
            (void)wpp.faceNormals();

            if (pp.size())
            {
                hasWedgePatches = true;
                wedgeDirVec += cmptMag(wpp.centreNormal());
            }
        }
    }


    if (returnReduceOr(hasEmptyPatches))
    {
        reduce(emptyDirVec, sumOp<vector>());

        emptyDirVec.normalise();

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (emptyDirVec[cmpt] > 1e-6)
            {
                solutionD_[cmpt] = -1;
            }
            else
            {
                solutionD_[cmpt] = 1;
            }
        }
    }


    // Knock out wedge directions

    geometricD_ = solutionD_;

    if (returnReduceOr(hasWedgePatches))
    {
        reduce(wedgeDirVec, sumOp<vector>());

        wedgeDirVec.normalise();

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (wedgeDirVec[cmpt] > 1e-6)
            {
                geometricD_[cmpt] = -1;
            }
            else
            {
                geometricD_[cmpt] = 1;
            }
        }
    }
}


Foam::autoPtr<Foam::labelIOList> Foam::polyMesh::readTetBasePtIs() const
{
    IOobject io
    (
        "tetBasePtIs",
        instance(),
        meshSubDir,
        *this,
        IOobject::LAZY_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (io.typeHeaderOk<labelIOList>(true))
    {
        return autoPtr<labelIOList>::New(io);
    }

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyMesh::polyMesh(const IOobject& io, const bool doInit)
:
    objectRegistry(io),
    primitiveMesh(),
    data_(static_cast<const objectRegistry&>(*this)),
    points_
    (
        IOobject
        (
            "points",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    faces_
    (
        IOobject
        (
            "faces",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    owner_
    (
        IOobject
        (
            "owner",
            faces_.instance(),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            faces_.instance(),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            time().findInstance // allow 'newer' boundary file
            (
                meshDir(),
                "boundary",
                IOobject::MUST_READ,
                faces_.instance()
            ),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    bounds_(points_),
    comm_(UPstream::worldComm),
    geometricD_(Zero),
    solutionD_(Zero),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            faces_.instance(),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        PtrList<entry>()
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            faces_.instance(),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        PtrList<entry>()
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            faces_.instance(),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this,
        PtrList<entry>()
    ),
    moving_(false),
    topoChanging_(false),
    storeOldCellCentres_(false),
    curMotionTimeIndex_(time().timeIndex())
{
    if (owner_.hasHeaderClass())
    {
        initMesh();
    }
    else
    {
        cellCompactIOList cLst
        (
            IOobject
            (
                "cells",
                time().findInstance(meshDir(), "cells"),
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Set the primitive mesh
        initMesh(cLst);

        owner_.write();
        neighbour_.write();
    }

    if (returnReduceOr(boundary_.empty()))
    {
        WarningInFunction
            << "Missing mesh boundary on one or more domains" << endl;

        // Warn if global empty mesh
        if (returnReduceAnd(!nPoints()))
        {
            WarningInFunction
                << "No points in mesh" << endl;
        }
        if (returnReduceAnd(!nCells()))
        {
            WarningInFunction
                << "No cells in mesh" << endl;
        }
    }

    if (doInit)
    {
        polyMesh::init(false);  // do not init lower levels
    }
}


bool Foam::polyMesh::init(const bool doInit)
{
    if (doInit)
    {
        primitiveMesh::init(doInit);
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    // Initialise demand-driven data
    calcDirections();

    return false;
}


Foam::polyMesh::polyMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& owner,
    labelList&& neighbour,
    const bool syncPar
)
:
    objectRegistry(io),
    primitiveMesh(),
    data_(static_cast<const objectRegistry&>(*this)),
    points_
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,  //io.readOpt(),
            io.writeOpt()
        ),
        std::move(points)
    ),
    faces_
    (
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,  //io.readOpt(),
            io.writeOpt()
        ),
        std::move(faces)
    ),
    owner_
    (
        IOobject
        (
            "owner",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,  //io.readOpt(),
            io.writeOpt()
        ),
        std::move(owner)
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,  //io.readOpt(),
            io.writeOpt()
        ),
        std::move(neighbour)
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,  // ignore since no alternative can be supplied
            io.writeOpt()
        ),
        *this,
        polyPatchList()
    ),
    bounds_(points_, syncPar),
    comm_(UPstream::worldComm),
    geometricD_(Zero),
    solutionD_(Zero),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,  // ignore since no alternative can be supplied
            IOobject::NO_WRITE
        ),
        *this,
        Foam::zero{}
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,  // ignore since no alternative can be supplied
            IOobject::NO_WRITE
        ),
        *this,
        Foam::zero{}
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,  // ignore since no alternative can be supplied
            IOobject::NO_WRITE
        ),
        *this,
        Foam::zero{}
    ),
    moving_(false),
    topoChanging_(false),
    storeOldCellCentres_(false),
    curMotionTimeIndex_(time().timeIndex())
{
    // Check if the faces and cells are valid
    forAll(faces_, facei)
    {
        const face& curFace = faces_[facei];

        if (min(curFace) < 0 || max(curFace) > points_.size())
        {
            FatalErrorInFunction
                << "Face " << facei << "contains vertex labels out of range: "
                << curFace << " Max point index = " << points_.size()
                << abort(FatalError);
        }
    }

    // Set the primitive mesh
    initMesh();
}


Foam::polyMesh::polyMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
:
    objectRegistry(io),
    primitiveMesh(),
    data_(static_cast<const objectRegistry&>(*this)),
    points_
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            io.writeOpt()
        ),
        std::move(points)
    ),
    faces_
    (
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            io.writeOpt()
        ),
        std::move(faces)
    ),
    owner_
    (
        IOobject
        (
            "owner",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            io.writeOpt()
        ),
        Foam::zero{}
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            io.writeOpt()
        ),
        Foam::zero{}
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            io.writeOpt()
        ),
        *this,
        Foam::zero{}
    ),
    bounds_(points_, syncPar),
    comm_(UPstream::worldComm),
    geometricD_(Zero),
    solutionD_(Zero),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        Foam::zero{}
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        Foam::zero{}
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        Foam::zero{}
    ),
    moving_(false),
    topoChanging_(false),
    storeOldCellCentres_(false),
    curMotionTimeIndex_(time().timeIndex())
{
    // Check if faces are valid
    forAll(faces_, facei)
    {
        const face& curFace = faces_[facei];

        if (min(curFace) < 0 || max(curFace) > points_.size())
        {
            FatalErrorInFunction
                << "Face " << facei << "contains vertex labels out of range: "
                << curFace << " Max point index = " << points_.size()
                << abort(FatalError);
        }
    }

    // Transfer in cell list
    cellList cLst(std::move(cells));

    // Check if cells are valid
    forAll(cLst, celli)
    {
        const cell& curCell = cLst[celli];

        if (min(curCell) < 0 || max(curCell) > faces_.size())
        {
            FatalErrorInFunction
                << "Cell " << celli << "contains face labels out of range: "
                << curCell << " Max face index = " << faces_.size()
                << abort(FatalError);
        }
    }

    // Set the primitive mesh
    initMesh(cLst);
}


Foam::polyMesh::polyMesh
(
    const IOobject& io,
    const Foam::zero,
    const bool syncPar
)
:
    polyMesh(io, pointField(), faceList(), labelList(), labelList(), syncPar)
{}


void Foam::polyMesh::resetPrimitives
(
    autoPtr<pointField>&& points,
    autoPtr<faceList>&& faces,
    autoPtr<labelList>&& owner,
    autoPtr<labelList>&& neighbour,
    const labelUList& patchSizes,
    const labelUList& patchStarts,
    const bool validBoundary
)
{
    // Clear addressing. Keep geometric props and updateable props for mapping.
    clearAddressing(true);

    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (points)
    {
        points_.transfer(*points);
        bounds_ = boundBox(points_, validBoundary);
    }

    if (faces)
    {
        faces_.transfer(*faces);
    }

    if (owner)
    {
        owner_.transfer(*owner);
    }

    if (neighbour)
    {
        neighbour_.transfer(*neighbour);
    }


    // Reset patch sizes and starts
    forAll(boundary_, patchi)
    {
        boundary_[patchi] = polyPatch
        (
            boundary_[patchi],
            boundary_,
            patchi,
            patchSizes[patchi],
            patchStarts[patchi]
        );
    }


    // Flags the mesh files as being changed
    setInstance(time().timeName());

    // Check if the faces and cells are valid
    forAll(faces_, facei)
    {
        const face& curFace = faces_[facei];

        if (min(curFace) < 0 || max(curFace) > points_.size())
        {
            FatalErrorInFunction
                << "Face " << facei << " contains vertex labels out of range: "
                << curFace << " Max point index = " << points_.size()
                << abort(FatalError);
        }
    }


    // Set the primitive mesh from the owner_, neighbour_.
    // Works out from patch end where the active faces stop.
    initMesh();


    if (validBoundary)
    {
        // Note that we assume that all the patches stay the same and are
        // correct etc. so we can already use the patches to do
        // processor-processor comms.

        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();

        // Warn if global empty mesh
        if (returnReduceAnd(!nPoints()) || returnReduceAnd(!nCells()))
        {
            FatalErrorInFunction
                << "No points or no cells in mesh" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyMesh::~polyMesh()
{
    clearOut();
    resetMotion();
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::word& Foam::polyMesh::regionName(const word& region)
{
    return (region == polyMesh::defaultRegion ? word::null : region);
}


Foam::fileName Foam::polyMesh::meshDir(const word& region)
{
    if (region.empty() || region == polyMesh::defaultRegion)
    {
        return polyMesh::meshSubDir;
    }

    return (region / polyMesh::meshSubDir);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fileName& Foam::polyMesh::dbDir() const
{
    if (objectRegistry::name() == polyMesh::defaultRegion)
    {
        return parent().dbDir();
    }

    return objectRegistry::dbDir();
}


Foam::fileName Foam::polyMesh::meshDir() const
{
    return dbDir()/meshSubDir;
}


const Foam::word& Foam::polyMesh::regionName() const
{
    return polyMesh::regionName(objectRegistry::name());
}


const Foam::fileName& Foam::polyMesh::pointsInstance() const
{
    return points_.instance();
}


const Foam::fileName& Foam::polyMesh::facesInstance() const
{
    return faces_.instance();
}


const Foam::Vector<Foam::label>& Foam::polyMesh::geometricD() const
{
    if (geometricD_.x() == 0)
    {
        calcDirections();
    }

    return geometricD_;
}


Foam::label Foam::polyMesh::nGeometricD() const
{
    return cmptSum(geometricD() + Vector<label>::one)/2;
}


const Foam::Vector<Foam::label>& Foam::polyMesh::solutionD() const
{
    if (solutionD_.x() == 0)
    {
        calcDirections();
    }

    return solutionD_;
}


Foam::label Foam::polyMesh::nSolutionD() const
{
    return cmptSum(solutionD() + Vector<label>::one)/2;
}


const Foam::labelIOList& Foam::polyMesh::tetBasePtIs() const
{
    if (!tetBasePtIsPtr_)
    {
        if (debug)
        {
            WarningInFunction
                << "Forcing storage of base points."
                << endl;
        }

        labelList basePts
        (
            polyMeshTetDecomposition::findFaceBasePts(*this)
        );

        tetBasePtIsPtr_.reset
        (
            new labelIOList
            (
                IOobject
                (
                    "tetBasePtIs",
                    instance(),
                    meshSubDir,
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                std::move(basePts)
            )
        );
    }

    return *tetBasePtIsPtr_;
}


const Foam::indexedOctree<Foam::treeDataCell>&
Foam::polyMesh::cellTree() const
{
    if (!cellTreePtr_)
    {
        Random rndGen(261782);

        treeBoundBox overallBb(points());
        overallBb.inflate(rndGen, 1e-4, ROOTVSMALL);

        cellTreePtr_.reset
        (
            new indexedOctree<treeDataCell>
            (
                treeDataCell
                (
                    false,      // not cache bb
                    *this,
                    CELL_TETS   // use tet-decomposition for any inside test
                ),
                overallBb,
                8,              // maxLevel
                10,             // leafsize
                5.0             // duplicity
            )
        );
    }

    return *cellTreePtr_;
}


void Foam::polyMesh::addPatches
(
    polyPatchList& plist,
    const bool validBoundary
)
{
    if (boundaryMesh().size())
    {
        FatalErrorInFunction
            << "boundary already exists"
            << abort(FatalError);
    }

    // Reset valid directions
    geometricD_ = Zero;
    solutionD_ = Zero;

    boundary_.transfer(plist);

    // parallelData depends on the processorPatch ordering so force
    // recalculation. Problem: should really be done in removeBoundary but
    // there is some info in parallelData which might be interesting inbetween
    // removeBoundary and addPatches.
    globalMeshDataPtr_.reset(nullptr);

    if (validBoundary)
    {
        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();

        boundary_.checkDefinition();
    }
}


void Foam::polyMesh::addZones
(
    PtrList<pointZone>&& pz,
    PtrList<faceZone>&& fz,
    PtrList<cellZone>&& cz
)
{
    if (pointZones_.size() || faceZones_.size() || cellZones_.size())
    {
        FatalErrorInFunction
            << "point, face or cell zone already exists"
            << abort(FatalError);
    }

    // Point zones - take ownership of the pointers
    if (pz.size())
    {
        pointZones_.clear();
        pointZones_.transfer(pz);
        pointZones_.writeOpt(IOobject::AUTO_WRITE);
    }

    // Face zones - take ownership of the pointers
    if (fz.size())
    {
        faceZones_.clear();
        faceZones_.transfer(fz);
        faceZones_.writeOpt(IOobject::AUTO_WRITE);
    }

    // Cell zones - take ownership of the pointers
    if (cz.size())
    {
        cellZones_.clear();
        cellZones_.transfer(cz);
        cellZones_.writeOpt(IOobject::AUTO_WRITE);
    }
}


void Foam::polyMesh::addPatches
(
    const List<polyPatch*>& p,
    const bool validBoundary
)
{
    // Acquire ownership of the pointers
    polyPatchList plist(const_cast<List<polyPatch*>&>(p));

    addPatches(plist, validBoundary);
}


void Foam::polyMesh::addZones
(
    const List<pointZone*>& pz,
    const List<faceZone*>& fz,
    const List<cellZone*>& cz
)
{
    // Acquire ownership of the pointers
    addZones
    (
        PtrList<pointZone>(const_cast<List<pointZone*>&>(pz)),
        PtrList<faceZone>(const_cast<List<faceZone*>&>(fz)),
        PtrList<cellZone>(const_cast<List<cellZone*>&>(cz))
    );
}


const Foam::pointField& Foam::polyMesh::points() const
{
    if (clearedPrimitives_)
    {
        FatalErrorInFunction
            << "points deallocated"
            << abort(FatalError);
    }

    return points_;
}


bool Foam::polyMesh::upToDatePoints(const regIOobject& io) const
{
    return io.upToDate(points_);
}


void Foam::polyMesh::setUpToDatePoints(regIOobject& io) const
{
    io.eventNo() = points_.eventNo()+1;
}


const Foam::faceList& Foam::polyMesh::faces() const
{
    if (clearedPrimitives_)
    {
        FatalErrorInFunction
            << "faces deallocated"
            << abort(FatalError);
    }

    return faces_;
}


const Foam::labelList& Foam::polyMesh::faceOwner() const
{
    return owner_;
}


const Foam::labelList& Foam::polyMesh::faceNeighbour() const
{
    return neighbour_;
}


const Foam::pointField& Foam::polyMesh::oldPoints() const
{
    if (!moving_)
    {
        return points_;
    }

    if (!oldPointsPtr_)
    {
        if (debug)
        {
            WarningInFunction << endl;
        }

        oldPointsPtr_.reset(new pointField(points_));
        curMotionTimeIndex_ = time().timeIndex();
    }

    return *oldPointsPtr_;
}


const Foam::pointField& Foam::polyMesh::oldCellCentres() const
{
    storeOldCellCentres_ = true;

    if (!moving_)
    {
        return cellCentres();
    }

    if (!oldCellCentresPtr_)
    {
        oldCellCentresPtr_.reset(new pointField(cellCentres()));
    }

    return *oldCellCentresPtr_;
}


void Foam::polyMesh::movePoints(const pointField& newPoints)
{
    DebugInFunction
        << "Moving points for time " << time().value()
        << " index " << time().timeIndex() << endl;

    if (newPoints.size() != points_.size())
    {
        FatalErrorInFunction
            << "Size of newPoints " << newPoints.size()
            << " does not correspond to current mesh points size "
            << points_.size()
            << exit(FatalError);
    }


    moving(true);

    // Pick up old points
    if (curMotionTimeIndex_ != time().timeIndex())
    {
        if (debug)
        {
            Info<< "void polyMesh::movePoints(const pointField&) : "
                << " Storing current points for time " << time().value()
                << " index " << time().timeIndex() << endl;
        }

        if (storeOldCellCentres_)
        {
            oldCellCentresPtr_.reset(nullptr);
            oldCellCentresPtr_.reset(new pointField(cellCentres()));
        }

        // Mesh motion in the new time step
        oldPointsPtr_.reset(nullptr);
        oldPointsPtr_.reset(new pointField(points_));
        curMotionTimeIndex_ = time().timeIndex();
    }

    points_ = newPoints;

    bool moveError = false;
    if (debug)
    {
        // Check mesh motion
        if (checkMeshMotion(points_, true))
        {
            moveError = true;

            InfoInFunction
                << "Moving the mesh with given points will "
                << "invalidate the mesh." << nl
                << "Mesh motion should not be executed." << endl;
        }
    }

    points_.writeOpt(IOobject::AUTO_WRITE);
    points_.instance() = time().timeName();
    points_.eventNo() = getEvent();

    if (tetBasePtIsPtr_)
    {
        tetBasePtIsPtr_->writeOpt(IOobject::AUTO_WRITE);
        tetBasePtIsPtr_->instance() = time().timeName();
        tetBasePtIsPtr_->eventNo() = getEvent();
    }

    // Currently a no-op; earlier versions set meshPhi and call
    // primitiveMesh::clearGeom
    (void)primitiveMesh::movePoints(points_, oldPoints());

    // Update the mesh geometry (via fvGeometryScheme)
    // - updateGeom is virtual -> calls fvMesh::updateGeom (or higher)
    // - fvMesh::updateGeom defers to surfaceInterpolation::updateGeom(),
    //   which defers to fvGeometryScheme::movePoints()
    // - set the mesh flux
    // - clear out/recalculate stale geometry
    updateGeom();

    // Adjust parallel shared points
    if (globalMeshDataPtr_)
    {
        globalMeshDataPtr_->movePoints(points_);
    }

    // Force recalculation of all geometric data with new points

    bounds_ = boundBox(points_);
    boundary_.movePoints(points_);

    pointZones_.movePoints(points_);
    faceZones_.movePoints(points_);
    cellZones_.movePoints(points_);

    // Reset cell tree - it gets built from mesh geometry so might have
    // wrong boxes. It is correct as long as none of the cells leaves
    // the boxes it is in which most likely is almost never the case except
    // for tiny displacements. An alternative is to check the displacements
    // to see if they are tiny - imagine a big windtunnel with a small rotating
    // object. In this case the processors without the rotating object wouldn't
    // have to clear any geometry. However your critical path still stays the
    // same so no time would be gained (unless the decomposition gets weighted).
    // Small benefit for lots of scope for problems so not done.
    cellTreePtr_.reset(nullptr);

    // Reset valid directions (could change with rotation)
    geometricD_ = Zero;
    solutionD_ = Zero;

    // Note: tet-base decomposition does not get cleared. Ideally your face
    // decomposition should not change during mesh motion ...


    meshObject::movePoints<polyMesh>(*this);
    meshObject::movePoints<pointMesh>(*this);

    const_cast<Time&>(time()).functionObjects().movePoints(*this);


    if (debug && moveError)
    {
        // Write mesh to ease debugging. Note we want to avoid calling
        // e.g. fvMesh::write since meshPhi not yet complete.
        polyMesh::write();
    }
}


void Foam::polyMesh::resetMotion() const
{
    curMotionTimeIndex_ = 0;
    oldPointsPtr_.reset(nullptr);
    oldCellCentresPtr_.reset(nullptr);
}


bool Foam::polyMesh::hasGlobalData() const noexcept
{
    return bool(globalMeshDataPtr_);
}


const Foam::globalMeshData& Foam::polyMesh::globalData() const
{
    if (!globalMeshDataPtr_)
    {
        if (debug)
        {
            Pout<< "polyMesh::globalData() const : "
                << "Constructing parallelData from processor topology"
                << endl;
        }
        // Construct globalMeshData using processorPatch information only.
        globalMeshDataPtr_.reset(new globalMeshData(*this));
    }

    return *globalMeshDataPtr_;
}


void Foam::polyMesh::removeFiles(const fileName& instanceDir) const
{
    fileName meshFilesPath = thisDb().time().path()/instanceDir/meshDir();

    rm(meshFilesPath/"points");
    rm(meshFilesPath/"faces");
    rm(meshFilesPath/"owner");
    rm(meshFilesPath/"neighbour");
    rm(meshFilesPath/"cells");
    rm(meshFilesPath/"boundary");
    rm(meshFilesPath/"pointZones");
    rm(meshFilesPath/"faceZones");
    rm(meshFilesPath/"cellZones");
    rm(meshFilesPath/"meshModifiers");
    rm(meshFilesPath/"parallelData");

    // remove subdirectories
    if (isDir(meshFilesPath/"sets"))
    {
        rmDir(meshFilesPath/"sets");
    }
}


void Foam::polyMesh::removeFiles() const
{
    removeFiles(instance());
}


void Foam::polyMesh::findCellFacePt
(
    const point& p,
    label& celli,
    label& tetFacei,
    label& tetPti
) const
{
    celli = -1;
    tetFacei = -1;
    tetPti = -1;

    const indexedOctree<treeDataCell>& tree = cellTree();

    // Find point inside cell
    celli = tree.findInside(p);

    if (celli != -1)
    {
        // Check the nearest cell to see if the point is inside.
        findTetFacePt(celli, p, tetFacei, tetPti);
    }
}


void Foam::polyMesh::findTetFacePt
(
    const label celli,
    const point& p,
    label& tetFacei,
    label& tetPti
) const
{
    const polyMesh& mesh = *this;

    tetIndices tet(polyMeshTetDecomposition::findTet(mesh, celli, p));
    tetFacei = tet.face();
    tetPti = tet.tetPt();
}


bool Foam::polyMesh::pointInCell
(
    const point& p,
    label celli,
    const cellDecomposition decompMode
) const
{
    switch (decompMode)
    {
        case FACE_PLANES:
        {
            return primitiveMesh::pointInCell(p, celli);
        }
        break;

        case FACE_CENTRE_TRIS:
        {
            // only test that point is on inside of plane defined by cell face
            // triangles
            const cell& cFaces = cells()[celli];

            forAll(cFaces, cFacei)
            {
                label facei = cFaces[cFacei];
                const face& f = faces_[facei];
                const point& fc = faceCentres()[facei];
                bool isOwn = (owner_[facei] == celli);

                forAll(f, fp)
                {
                    label pointi;
                    label nextPointi;

                    if (isOwn)
                    {
                        pointi = f[fp];
                        nextPointi = f.nextLabel(fp);
                    }
                    else
                    {
                        pointi = f.nextLabel(fp);
                        nextPointi = f[fp];
                    }

                    triPointRef faceTri
                    (
                        points()[pointi],
                        points()[nextPointi],
                        fc
                    );

                    vector proj = p - faceTri.centre();

                    if ((faceTri.areaNormal() & proj) > 0)
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        break;

        case FACE_DIAG_TRIS:
        {
            // only test that point is on inside of plane defined by cell face
            // triangles
            const cell& cFaces = cells()[celli];

            forAll(cFaces, cFacei)
            {
                label facei = cFaces[cFacei];
                const face& f = faces_[facei];

                for (label tetPti = 1; tetPti < f.size() - 1; tetPti++)
                {
                    // Get tetIndices of face triangle
                    tetIndices faceTetIs(celli, facei, tetPti);

                    triPointRef faceTri = faceTetIs.faceTri(*this);

                    vector proj = p - faceTri.centre();

                    if ((faceTri.areaNormal() & proj) > 0)
                    {
                        return false;
                    }
                }
            }

            return true;
        }
        break;

        case CELL_TETS:
        {
            label tetFacei;
            label tetPti;

            findTetFacePt(celli, p, tetFacei, tetPti);

            return tetFacei != -1;
        }
        break;
    }

    return false;
}


Foam::label Foam::polyMesh::findCell
(
    const point& p,
    const cellDecomposition decompMode
) const
{
    if
    (
        Pstream::parRun()
     && (decompMode == FACE_DIAG_TRIS || decompMode == CELL_TETS)
    )
    {
        // Force construction of face-diagonal decomposition before testing
        // for zero cells.
        //
        // If parallel running a local domain might have zero cells so never
        // construct the face-diagonal decomposition which uses parallel
        // transfers.
        (void)tetBasePtIs();
    }

    if (nCells() == 0)
    {
        return -1;
    }

    if (decompMode == CELL_TETS)
    {
        // Advanced search method utilizing an octree
        // and tet-decomposition of the cells

        label celli;
        label tetFacei;
        label tetPti;

        findCellFacePt(p, celli, tetFacei, tetPti);

        return celli;
    }
    else
    {
        // Approximate search avoiding the construction of an octree
        // and cell decomposition

        if (Pstream::parRun() && decompMode == FACE_DIAG_TRIS)
        {
            // Force construction of face-diagonal decomposition before testing
            // for zero cells. If parallel running a local domain might have
            // zero cells so never construct the face-diagonal decomposition
            // (which uses parallel transfers)
            (void)tetBasePtIs();
        }

        // Find the nearest cell centre to this location
        label celli = findNearestCell(p);

        // If point is in the nearest cell return
        if (pointInCell(p, celli, decompMode))
        {
            return celli;
        }
        else
        {
            // Point is not in the nearest cell so search all cells

            for (label celli = 0; celli < nCells(); celli++)
            {
                if (pointInCell(p, celli, decompMode))
                {
                    return celli;
                }
            }

            return -1;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::polyMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    // Currently no special treatment. Just write the objects

    return objectRegistry::writeObject(streamOpt, writeOnProc);
}


// ************************************************************************* //
