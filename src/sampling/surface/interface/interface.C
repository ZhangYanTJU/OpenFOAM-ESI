/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 DLR
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

#include "interface.H"
#include "reconstructionSchemes.H"
#include "cutCellPLIC.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interface, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interface::interface
(
    const fvMesh& mesh
)
:
    MeshStorage(),
    mesh_(mesh)
{
    reconstructionSchemes& surf =
        mesh_.lookupObjectRef<reconstructionSchemes>("reconstructionScheme");

    surf.reconstruct();

    cutCellPLIC cellCut(mesh_);

    const volVectorField& normal = surf.normal();
    const volVectorField& centre = surf.centre();

    const boolList& interfaceCells = surf.interfaceCell();

    DynamicList<List<point>> facePts;
    DynamicList<label> interfaceCellAdressing(0.1*mesh_.nCells());

    forAll(interfaceCells,cellI)
    {
        if (interfaceCells[cellI])
        {
            if (mag(normal[cellI]) != 0)
            {
                interfaceCellAdressing.append(cellI);
                vector n = -normal[cellI]/mag(normal[cellI]);

                scalar cutVal = (centre[cellI]-mesh_.C()[cellI]) & n;

                cellCut.calcSubCell(cellI,cutVal,n);
                facePts.append(cellCut.facePoints());
            }
        }
    }

    meshCells_.setSize(interfaceCellAdressing.size());

    forAll(meshCells_,i)
    {
        meshCells_[i] = interfaceCellAdressing[i];
    }

    // Transfer to mesh storage
    {
        faceList faces(facePts.size());

        label nPoints = 0;
        forAll(facePts,i)
        {
            face f(facePts[i].size());
            forAll(f,fi)
            {
                f[fi] = nPoints + fi;
            }
            faces[i] = f;

            nPoints += facePts[i].size();
        }
        pointField pts(nPoints);

        nPoints = 0; // reuse
        forAll(facePts,i)
        {
            forAll(facePts[i],fi)
            {
                pts[nPoints] = facePts[i][fi];
                ++nPoints;
            }
        }

        MeshStorage updated(std::move(pts), std::move(faces), surfZoneList());

        this->MeshStorage::transfer(updated);
    }
}


// ************************************************************************* //
