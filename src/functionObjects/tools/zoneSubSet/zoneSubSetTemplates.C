/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "volFieldsFwd.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::zoneSubSet::mapToZone
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volField;

    const labelList& cellMap = subSetMeshPtr_->cellMap();

    auto tsubSet = tmp<volField>::New
    (
        IOobject
        (
            "tsubSet",
            vf.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensioned<Type>(vf.dimensions())
    );
    auto& subSet = tsubSet.ref();

    // Map from sub-mesh to global mesh
    forAll (cellMap, celli)
    {
        subSet[cellMap[celli]] = vf[celli];
    }

    auto tcellZonesField = tmp<volField>::New
    (
        IOobject
        (
            vf.name(),
            vf.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensioned<Type>(vf.dimensions())
    );
    auto& cellZonesField = tcellZonesField.ref();

    // Map field in global mesh on original zones
    for (const word& zoneName : zoneNames_)
    {
        const labelList& cells = mesh_.cellZones()[zoneName];

        for (const label cell : cells)
        {
            cellZonesField[cell] = subSet[cell];
        }
    }

    return tcellZonesField;
}


// ************************************************************************* //
