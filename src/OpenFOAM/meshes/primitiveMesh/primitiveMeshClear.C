/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::primitiveMesh::printAllocated() const
{
    Pout<< "primitiveMesh allocated :" << endl;

    // Topology
    if (cellShapesPtr_)
    {
        Pout<< "    Cell shapes" << endl;
    }

    if (edgesPtr_)
    {
        Pout<< "    Edges" << endl;
    }

    if (ccPtr_)
    {
        Pout<< "    Cell-cells" << endl;
    }

    if (ecPtr_)
    {
        Pout<< "    Edge-cells" << endl;
    }

    if (pcPtr_)
    {
        Pout<< "    Point-cells" << endl;
    }

    if (cfPtr_)
    {
        Pout<< "    Cell-faces" << endl;
    }

    if (efPtr_)
    {
        Pout<< "    Edge-faces" << endl;
    }

    if (pfPtr_)
    {
        Pout<< "    Point-faces" << endl;
    }

    if (cePtr_)
    {
        Pout<< "    Cell-edges" << endl;
    }

    if (fePtr_)
    {
        Pout<< "    Face-edges" << endl;
    }

    if (pePtr_)
    {
        Pout<< "    Point-edges" << endl;
    }

    if (ppPtr_)
    {
        Pout<< "    Point-point" << endl;
    }

    if (cpPtr_)
    {
        Pout<< "    Cell-point" << endl;
    }

    // Geometry
    if (cellCentresPtr_)
    {
        Pout<< "    Cell-centres" << endl;
    }

    if (cellVolumesPtr_)
    {
        Pout<< "    Cell-volumes" << endl;
    }

    if (faceCentresPtr_)
    {
        Pout<< "    Face-centres" << endl;
    }

    if (faceAreasPtr_)
    {
        Pout<< "    Face-areas" << endl;
    }
}


void Foam::primitiveMesh::clearGeom()
{
    if (debug)
    {
        Pout<< "primitiveMesh::clearGeom() : "
            << "clearing geometric data"
            << endl;
    }

    cellCentresPtr_.reset(nullptr);
    cellVolumesPtr_.reset(nullptr);
    faceCentresPtr_.reset(nullptr);
    faceAreasPtr_.reset(nullptr);
}


void Foam::primitiveMesh::clearCellGeom()
{
    if (debug)
    {
        Pout<< "primitiveMesh::clearCellGeom() : "
            << "clearing cell centres and volumes"
            << endl;
    }

    cellCentresPtr_.reset(nullptr);
    cellVolumesPtr_.reset(nullptr);
}


void Foam::primitiveMesh::clearAddressing()
{
    if (debug)
    {
        Pout<< "primitiveMesh::clearAddressing() : "
            << "clearing topology"
            << endl;
    }

    cellShapesPtr_.reset(nullptr);

    clearOutEdges();

    ccPtr_.reset(nullptr);
    ecPtr_.reset(nullptr);
    pcPtr_.reset(nullptr);

    cfPtr_.reset(nullptr);
    efPtr_.reset(nullptr);
    pfPtr_.reset(nullptr);

    cePtr_.reset(nullptr);
    fePtr_.reset(nullptr);
    pePtr_.reset(nullptr);
    ppPtr_.reset(nullptr);
    cpPtr_.reset(nullptr);
}


void Foam::primitiveMesh::clearOut()
{
    clearGeom();
    clearAddressing();
}


// ************************************************************************* //
