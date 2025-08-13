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

#include "internalPointProber.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(internalPointProber, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::internalPointProber::findElements(const fvMesh& mesh)
{
    DebugInfo<< "internalPointProber: resetting sample locations" << endl;

    elementList_.resize_nocopy(pointField::size());
    faceList_.resize_nocopy(pointField::size());
    processor_.resize_nocopy(pointField::size());
    processor_ = -1;

    forAll(*this, probei)
    {
        const point& location = (*this)[probei];

        const label celli = mesh.findCell(location);

        elementList_[probei] = celli;

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
            faceList_[probei] = minFaceID;
        }
        else
        {
            faceList_[probei] = -1;
        }

        if (debug && (elementList_[probei] != -1 || faceList_[probei] != -1))
        {
            Pout<< "internalPointProber : found point " << location
                << " in cell " << elementList_[probei]
                << " and face " << faceList_[probei] << endl;
        }
    }


    // Check if all probes have been found.
    forAll(elementList_, probei)
    {
        const point& location = operator[](probei);
        label celli = elementList_[probei];
        label facei = faceList_[probei];

        processor_[probei] = (celli != -1 ? Pstream::myProcNo() : -1);

        // Check at least one processor with cell.
        reduce(celli, maxOp<label>());
        reduce(facei, maxOp<label>());
        reduce(processor_[probei], maxOp<label>());

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
            if (elementList_[probei] != -1 && elementList_[probei] != celli)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << elementList_[probei]
                    << " on my domain " << Pstream::myProcNo()
                    << " and cell " << celli << " on some other domain."
                    << nl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }

            if (faceList_[probei] != -1 && faceList_[probei] != facei)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << faceList_[probei]
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

Foam::internalPointProber::internalPointProber
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pointProberBase(mesh, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::internalPointProber::read(const dictionary& dict)
{
    if (!pointProberBase::read(dict))
    {
        return false;
    }

    // Initialise cells to sample from supplied locations
    findElements(mesh_);

    return true;
}


// ************************************************************************* //
