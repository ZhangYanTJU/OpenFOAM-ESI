/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "fvcCellReduce.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CombineOp>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvc::cellReduce
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const CombineOp& cop,
    const Type& nullValue
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    const fvMesh& mesh = ssf.mesh();

    tmp<volFieldType> tresult
    (
        new volFieldType
        (
            IOobject
            (
                "cellReduce(" + ssf.name() + ')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("initialValue", ssf.dimensions(), nullValue),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );

    volFieldType& result = tresult.ref();

    const labelUList& own = mesh.owner();
    const labelUList& nbr = mesh.neighbour();

    forAll(own, i)
    {
        const label celli = own[i];
        cop(result[celli], ssf[i]);
    }
    forAll(nbr, i)
    {
        const label celli = nbr[i];
        cop(result[celli], ssf[i]);
    }

    forAll(mesh.boundary(), patchi)
    {
        const auto& pFaceCells = mesh.boundary()[patchi].faceCells();
        const auto& pssf = ssf.boundaryField()[patchi];

        forAll(pssf, i)
        {
            const label celli = pFaceCells[i];
            cop(result[celli], pssf[i]);
        }
    }

    result.correctBoundaryConditions();

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvc::cellReduce
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf,
    const CombineOp& cop,
    const Type& nullValue
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        cellReduce
        (
            tssf,
            cop,
            nullValue
        )
    );

    tssf.clear();

    return tvf;
}


// ************************************************************************* //
