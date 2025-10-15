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

#include "prober.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(prober, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prober::prober
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    samplePointScheme_("cell")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::prober::read(const dictionary& dict)
{
    dict.readEntry("probeLocations", probes_);

    if (probes_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Empty 'probeLocations' list."
            << exit(FatalIOError);
    }

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


void Foam::prober::updateMesh(const mapPolyMesh& mpm)
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
            DynamicList<label> elems(cellIds_.size());

            const labelList& reverseMap = mpm.reverseCellMap();
            forAll(cellIds_, i)
            {
                label celli = cellIds_[i];
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

            cellIds_.transfer(elems);
        }

        // 2. Update faces
        {
            DynamicList<label> elems(faceIds_.size());

            const labelList& reverseMap = mpm.reverseFaceMap();
            for (const label facei : faceIds_)
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

            faceIds_.transfer(elems);
        }
    }
}


void Foam::prober::movePoints(const polyMesh& mesh)
{
    DebugInfo<< "probes: movePoints" << endl;

    if (fixedLocations_ && &mesh == &mesh_)
    {
        this->findElements(mesh_);
    }
}


// ************************************************************************* //
