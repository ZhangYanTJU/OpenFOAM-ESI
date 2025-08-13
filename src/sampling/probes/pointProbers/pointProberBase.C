/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "pointProberBase.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointProberBase, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointProberBase::pointProberBase
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pointField(),
    mesh_(mesh),
    samplePointScheme_("cell")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointProberBase::read(const dictionary& dict)
{
    dict.readEntry("probeLocations", static_cast<pointField&>(*this));

    fixedLocations_ = dict.getOrDefault<bool>("fixedLocations", true);
    includeOutOfBounds_ = dict.getOrDefault<bool>("includeOutOfBounds", true);

    if (dict.readIfPresent("interpolationScheme", samplePointScheme_))
    {
        if (!fixedLocations_ && samplePointScheme_ != "cell")
        {
            WarningInFunction
                << "Only cell interpolation can be applied when "
                << "not using fixedLocations. InterpolationScheme "
                << "entry will be ignored"
                << endl;
        }
    }

    return true;
}


void Foam::pointProberBase::updateMesh(const mapPolyMesh& mpm)
{
    DebugInfo<< "probes: updateMesh" << endl;

    if (&mpm.mesh() != &mesh_)
    {
        return;
    }

    if (fixedLocations_)
    {
        this->findElements(mesh_);
    }
    else
    {
        DebugInfo<< "probes: remapping sample locations" << endl;

        // 1. Update cells
        {
            DynamicList<label> elems(elementList_.size());

            const labelList& reverseMap = mpm.reverseCellMap();
            forAll(elementList_, i)
            {
                label celli = elementList_[i];
                if (celli != -1)
                {
                    label newCelli = reverseMap[celli];
                    if (newCelli == -1)
                    {
                        // cell removed
                    }
                    else if (newCelli < -1)
                    {
                        // cell merged
                        elems.append(-newCelli - 2);
                    }
                    else
                    {
                        // valid new cell
                        elems.append(newCelli);
                    }
                }
                else
                {
                    // Keep -1 elements so the size stays the same
                    elems.append(-1);
                }
            }

            elementList_.transfer(elems);
        }

        // 2. Update faces
        {
            DynamicList<label> elems(faceList_.size());

            const labelList& reverseMap = mpm.reverseFaceMap();
            for (const label facei : faceList_)
            {
                if (facei != -1)
                {
                    label newFacei = reverseMap[facei];
                    if (newFacei == -1)
                    {
                        // face removed
                    }
                    else if (newFacei < -1)
                    {
                        // face merged
                        elems.append(-newFacei - 2);
                    }
                    else
                    {
                        // valid new face
                        elems.append(newFacei);
                    }
                }
                else
                {
                    // Keep -1 elements
                    elems.append(-1);
                }
            }

            faceList_.transfer(elems);
        }
    }
}


void Foam::pointProberBase::movePoints(const polyMesh& mesh)
{
    DebugInfo<< "probes: movePoints" << endl;

    if (fixedLocations_ && &mesh == &mesh_)
    {
        this->findElements(mesh_);
    }
}


// ************************************************************************* //
