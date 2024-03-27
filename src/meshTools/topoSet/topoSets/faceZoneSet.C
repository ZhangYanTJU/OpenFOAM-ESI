/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "faceZoneSet.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "setToFaceZone.H"
#include "setsToFaceZone.H"
#include "syncTools.H"
#include "ListOps.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZoneSet, 0);
    addToRunTimeSelectionTable(topoSet, faceZoneSet, word);
    addToRunTimeSelectionTable(topoSet, faceZoneSet, size);
    addToRunTimeSelectionTable(topoSet, faceZoneSet, set);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceZoneSet::updateSet()
{
    labelList order(sortedOrder(addressing_));
    addressing_ = labelUIndList(addressing_, order)();
    flipMap_ = boolUIndList(flipMap_, order)();

    faceSet::clearStorage();
    faceSet::reserve(addressing_.size());
    faceSet::set(addressing_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZoneSet::faceZoneSet
(
    const polyMesh& mesh,
    const word& name,
    IOobjectOption::readOption rOpt,
    IOobjectOption::writeOption wOpt
)
:
    faceSet(mesh, name, 1024),  // do not read faceSet
    mesh_(mesh),
    addressing_(),
    flipMap_()
{
    const faceZoneMesh& faceZones = mesh.faceZones();
    label zoneID = faceZones.findZoneID(name);

    if (IOobjectOption::isReadRequired(rOpt) && zoneID == -1)
    {
        FatalErrorInFunction
            << "Zone named " << name << " not found.  "
            << "List of available zone names: " << faceZones.names()
            << exit(FatalError);
    }

    if
    (
         IOobjectOption::isReadRequired(rOpt)
     || (IOobjectOption::isReadOptional(rOpt) && zoneID != -1)
    )
    {
        const faceZone& fz = faceZones[zoneID];
        addressing_ = fz.addressing();
        flipMap_ = fz.flipMap();
    }

    updateSet();

    check(mesh.nFaces());
}


Foam::faceZoneSet::faceZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    IOobjectOption::writeOption wOpt
)
:
    faceSet(mesh, name, size, wOpt),
    mesh_(mesh),
    addressing_(),
    flipMap_()
{
    updateSet();
}


Foam::faceZoneSet::faceZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    IOobjectOption::writeOption wOpt
)
:
    faceSet(mesh, name, set.size(), wOpt),
    mesh_(mesh),
    addressing_(refCast<const faceZoneSet>(set).addressing()),
    flipMap_(refCast<const faceZoneSet>(set).flipMap())
{
    updateSet();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceZoneSet::invert(const label maxLen)
{
    // Count
    label n = 0;

    for (label facei = 0; facei < maxLen; ++facei)
    {
        if (!found(facei))
        {
            ++n;
        }
    }

    // Fill
    addressing_.setSize(n);
    flipMap_.setSize(n);
    n = 0;

    for (label facei = 0; facei < maxLen; ++facei)
    {
        if (!found(facei))
        {
            addressing_[n] = facei;
            flipMap_[n] = false;         //? or true?
            ++n;
        }
    }
    updateSet();
}


void Foam::faceZoneSet::subset
(
    const word& setName,
    const labelUList& setAddressing,
    const UList<bool>& setFlipMap
)
{
    label nConflict = 0;

    DynamicList<label> newAddressing(addressing_.size());
    DynamicList<bool> newFlipMap(flipMap_.size());

    Map<label> faceToIndex(invertToMap(addressing_));

    forAll(setAddressing, i)
    {
        const label facei = setAddressing[i];

        const auto iter = faceToIndex.cfind(facei);

        if (iter.good())
        {
            const label index = iter.val();

            if (setFlipMap.size() && (setFlipMap[i] != flipMap_[index]))
            {
                ++nConflict;
            }
            newAddressing.append(facei);
            newFlipMap.append(flipMap_[index]);
        }
    }

    if (nConflict > 0)
    {
        WarningInFunction
            << "subset : there are " << nConflict
            << " faces with different orientation in faceZoneSets "
            << name() << " and " << setName << endl;
    }

    addressing_.transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
    updateSet();
}


void Foam::faceZoneSet::subset(const topoSet& set)
{
    const auto* setPtr = dynamic_cast<const faceZoneSet*>(&set);

    if (setPtr)
    {
        subset(setPtr->name(), setPtr->addressing(), setPtr->flipMap());
    }
    else
    {
        // Assume a faceSet. Ignore flipMap
        const auto& fSet = refCast<const faceSet>(set);
        subset(fSet.name(), fSet.sortedToc(), boolList::null());
    }
}


void Foam::faceZoneSet::subset(const labelUList& set)
{
    subset(word::null, set, boolList::null());
}


void Foam::faceZoneSet::addSet
(
    const word& setName,
    const labelUList& setAddressing,
    const UList<bool>& setFlipMap
)
{
    label nConflict = 0;

    DynamicList<label> newAddressing(addressing_);
    DynamicList<bool> newFlipMap(flipMap_);

    Map<label> faceToIndex(invertToMap(addressing_));

    forAll(setAddressing, i)
    {
        const label facei = setAddressing[i];
        const auto iter = faceToIndex.cfind(facei);

        if (iter.good())
        {
            const label index = iter.val();

            if (setFlipMap.size() && (setFlipMap[i] != flipMap_[index]))
            {
                ++nConflict;
            }
        }
        else
        {
            newAddressing.append(facei);
            newFlipMap.append(setFlipMap.size() ? setFlipMap[i] : false);
        }
    }

    if (nConflict > 0)
    {
        WarningInFunction
            << "addSet : there are " << nConflict
            << " faces with different orientation in faceZonesSets "
            << name() << " and " << setName << endl;
    }

    addressing_.transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
    updateSet();
}


void Foam::faceZoneSet::addSet(const topoSet& set)
{
    const auto* setPtr = dynamic_cast<const faceZoneSet*>(&set);

    if (setPtr)
    {
        addSet(setPtr->name(), setPtr->addressing(), setPtr->flipMap());
    }
    else
    {
        // Assume a faceSet. Ignore flipMap
        const auto& fSet = refCast<const faceSet>(set);
        addSet(fSet.name(), fSet.sortedToc(), boolList::null());
    }
}


void Foam::faceZoneSet::addSet(const labelUList& set)
{
    addSet(word::null, set, boolList::null());
}


void Foam::faceZoneSet::subtractSet
(
    const word& setName,
    const labelUList& setAddressing,
    const UList<bool>& setFlipMap
)
{
    label nConflict = 0;

    DynamicList<label> newAddressing(addressing_.size());
    DynamicList<bool> newFlipMap(flipMap_.size());

    Map<label> faceToIndex(invertToMap(setAddressing));

    forAll(addressing_, i)
    {
        const label facei = addressing_[i];

        const auto iter = faceToIndex.cfind(facei);

        if (iter.good())
        {
            const label index = iter.val();

            if (setFlipMap.size() && (setFlipMap[index] != flipMap_[i]))
            {
                ++nConflict;
            }
        }
        else
        {
            // Not found in zoneSet so add
            newAddressing.append(facei);
            newFlipMap.append(setFlipMap.size() ? setFlipMap[i] : false);
        }
    }

    if (nConflict > 0)
    {
        WarningInFunction
            << "subtractSet : there are " << nConflict
            << " faces with different orientation in faceZonesSets "
            << name() << " and " << setName << endl;
    }

    addressing_.transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
    updateSet();
}


void Foam::faceZoneSet::subtractSet(const topoSet& set)
{
    const auto* setPtr = dynamic_cast<const faceZoneSet*>(&set);

    if (setPtr)
    {
        subtractSet(setPtr->name(), setPtr->addressing(), setPtr->flipMap());
    }
    else
    {
        // Assume a faceSet. Ignore flipMap
        const auto& fSet = refCast<const faceSet>(set);
        subtractSet(fSet.name(), fSet.sortedToc(), boolList::null());
    }
}


void Foam::faceZoneSet::subtractSet(const labelUList& set)
{
    subtractSet(word::null, set, boolList::null());
}


void Foam::faceZoneSet::sync(const polyMesh& mesh)
{
    // This routine serves two purposes
    // 1. make sure that any previous faceZoneSet manipulation is
    //    consistent across coupled boundaries
    // 2. push faceZone contents to faceSet (looses flip bit)


    // Collect all current zone info
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 0 : not in faceZone
    // 1 : in faceZone and unflipped
    //-1 : in faceZone and flipped
    const label UNFLIPPED = 1;
    const label FLIPPED = -1;
    labelList myZoneFace(mesh.nFaces(), Zero);

    forAll(addressing_, i)
    {
        const label facei = addressing_[i];
        myZoneFace[facei] =
        (
            flipMap_[i]
          ? FLIPPED
          : UNFLIPPED
        );
    }

    labelList neiZoneFace
    (
        SubList<label>
        (
            myZoneFace,
            mesh.nBoundaryFaces(),
            mesh.nInternalFaces()
        )
    );
    syncTools::swapBoundaryFaceList(mesh, neiZoneFace);


    const bitSet isMasterFace(syncTools::getMasterFaces(mesh));


    // Rebuild faceZone addressing and flipMap
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const labelHashSet& set = *this;

    DynamicList<label> newAddressing(set.size());
    DynamicList<bool> newFlipMap(set.size());

    for (const label facei : set)
    {
        // See if any info from original. If so maintain flipMap.
        if (facei < mesh.nInternalFaces())
        {
            newAddressing.append(facei);
            newFlipMap.append(myZoneFace[facei] == FLIPPED);
        }
        else
        {
            const label myStat = myZoneFace[facei];
            const label neiStat = neiZoneFace[facei-mesh.nInternalFaces()];

            if (myStat == 0)
            {
                // My face was not in zone. Check neighbour

                if (neiStat == UNFLIPPED)
                {
                    // Neighbour is unflipped so I am flipped
                    newAddressing.append(facei);
                    newFlipMap.append(true);
                }
                else if (neiStat == FLIPPED)
                {
                    newAddressing.append(facei);
                    newFlipMap.append(false);
                }
                else //if (neiStat == 0)
                {
                    // neighbour face not in zone either. Masterface decides.
                    newAddressing.append(facei);
                    newFlipMap.append(!isMasterFace[facei]);
                }
            }
            else
            {
                if (myStat == neiStat)
                {
                    // Conflict. masterFace wins
                    newAddressing.append(facei);
                    if (isMasterFace[facei])
                    {
                        newFlipMap.append(myStat == FLIPPED);
                    }
                    else
                    {
                        newFlipMap.append(neiStat == UNFLIPPED);
                    }
                }
                else
                {
                    newAddressing.append(facei);
                    newFlipMap.append(myStat == FLIPPED);
                }
            }
        }
    }

    addressing_.transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
    updateSet();
}


Foam::label Foam::faceZoneSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nFaces();
}


bool Foam::faceZoneSet::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    // Write shadow faceSet
    word oldTypeName = typeName;
    const_cast<word&>(type()) = faceSet::typeName;
    bool ok = faceSet::writeObject(streamOpt, writeOnProc);
    const_cast<word&>(type()) = oldTypeName;

    // Modify faceZone
    faceZoneMesh& faceZones = const_cast<polyMesh&>(mesh_).faceZones();
    label zoneID = faceZones.findZoneID(name());

    if (zoneID == -1)
    {
        zoneID = faceZones.size();

        faceZones.emplace_back
        (
            name(),
            addressing_,
            flipMap_,
            zoneID,
            faceZones
        );
    }
    else
    {
        faceZones[zoneID].resetAddressing(addressing_, flipMap_);
    }
    faceZones.clearAddressing();

    return ok && faceZones.write(writeOnProc);
}


void Foam::faceZoneSet::updateMesh(const mapPolyMesh& morphMap)
{
    // faceZone
    labelList newAddressing(addressing_.size());
    boolList newFlipMap(flipMap_.size(), false);

    label n = 0;
    forAll(addressing_, i)
    {
        label facei = addressing_[i];
        label newFacei = morphMap.reverseFaceMap()[facei];
        if (newFacei >= 0)
        {
            newAddressing[n] = newFacei;
            newFlipMap[n] = flipMap_[i];
            n++;
        }
    }
    newAddressing.setSize(n);
    newFlipMap.setSize(n);

    addressing_.transfer(newAddressing);
    flipMap_.transfer(newFlipMap);

    updateSet();
}


void Foam::faceZoneSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    faceSet::writeDebug(os, mesh, maxLen);
}


// ************************************************************************* //
