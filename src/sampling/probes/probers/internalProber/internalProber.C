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

#include "internalProber.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(internalProber, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::internalProber::findElements(const fvMesh& mesh)
{
    DebugInfo<< "internalProber: resetting sample locations" << endl;

    const pointField& probeLocations = this->probeLocations();

    cellIds_.resize_nocopy(probeLocations.size());
    faceIds_.resize_nocopy(probeLocations.size());
    procIds_.resize_nocopy(probeLocations.size());
    procIds_ = -1;

    forAll(probeLocations, probei)
    {
        const point& location = probeLocations[probei];

        const label celli = mesh.findCell(location);

        cellIds_[probei] = celli;

        if (celli != -1)
        {
            const labelList& cellFaces = mesh.cells()[celli];
            const vector& cellCentre = mesh.cellCentres()[celli];
            scalar minDistance = GREAT;
            label minFaceID = -1;
            forAll(cellFaces, i)
            {
                label facei = cellFaces[i];
                vector dist = mesh.faceCentres()[facei] - cellCentre;
                if (mag(dist) < minDistance)
                {
                    minDistance = mag(dist);
                    minFaceID = facei;
                }
            }
            faceIds_[probei] = minFaceID;
        }
        else
        {
            faceIds_[probei] = -1;
        }

        if (debug && (cellIds_[probei] != -1 || faceIds_[probei] != -1))
        {
            Pout<< "internalProber : found point " << location
                << " in cell " << cellIds_[probei]
                << " and face " << faceIds_[probei] << endl;
        }
    }


    // Check if all probes have been found.
    forAll(cellIds_, probei)
    {
        const point& location = probeLocations[probei];
        label celli = cellIds_[probei];
        label facei = faceIds_[probei];

        procIds_[probei] = (celli != -1 ? Pstream::myProcNo() : -1);

        // Check at least one processor with cell.
        reduce(celli, maxOp<label>());
        reduce(facei, maxOp<label>());
        reduce(procIds_[probei], maxOp<label>());

        if (celli == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any cell. Skipping location." << endl;
            }
        }
        else if (facei == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any face. Skipping location." << endl;
            }
        }
        else
        {
            // Make sure location not on two domains.
            if (cellIds_[probei] != -1 && cellIds_[probei] != celli)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << cellIds_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and cell " << celli << " on some other domain."
                    << nl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }

            if (faceIds_[probei] != -1 && faceIds_[probei] != facei)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << faceIds_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and face " << facei << " on some other domain."
                    << nl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::internalProber::internalProber
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    prober(mesh, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::internalProber::read(const dictionary& dict)
{
    if (!prober::read(dict))
    {
        return false;
    }

    // Initialise cells to sample from supplied locations
    findElements(mesh_);

    return true;
}


// ************************************************************************* //
