/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "assemblyDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "globalIndex.H"
#include "regionProperties.H"
#include "fvMesh.H"
#include "lduPrimitiveMeshAssembly.H"
#include "labelIOList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(assemblyDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        assemblyDecomp,
        dictionary
    );
}

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static inline fvMesh* newFvMesh(const word& meshName, const polyMesh& mesh)
{
    return new fvMesh
    (
        IOobject
        (
            meshName,
            mesh.time().timeName(),
            mesh.time(),
            IOobject::MUST_READ
        )
    );
}

} // End namespace fvMesh


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::assemblyDecomp::assemblyDecomp(const dictionary& decompDict)
:
    decompositionMethod(decompDict),
    methodDict_(findCoeffsDict(typeName + "Coeffs", selectionType::MANDATORY)),
    dataFile_
    (
        findCoeffsDict(typeName + "Coeffs").get<fileName>("dataFile")
    )
{
    method_ = decompositionMethod::New(methodDict_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::assemblyDecomp::decompose
(
    const polyMesh& mesh,
    const pointField&,
    const scalarField&
) const
{
    regionProperties rp(mesh.time());

    const wordList& fluidNames(rp["fluid"]);
    const wordList& solidsNames(rp["solid"]);

    // Use regionProperties to find all other regions
    UPtrList<fvMesh> meshes(fluidNames.size() + solidsNames.size());

    label regioni = 0;
    for (const word& regName : fluidNames)
    {
        meshes.set(regioni, newFvMesh(regName, mesh));
        ++regioni;
    }

    for (const word& regName : solidsNames)
    {
        meshes.set(regioni, newFvMesh(regName, mesh));
        ++regioni;
    }

    lduPrimitiveMeshAssembly
        assemblyLduMesh
        (
            meshes,
            IOobject
            (
                "assemblyLdu",
                mesh.time().timeName(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    const lduAddressing& addr = assemblyLduMesh.lduAddr();

    globalIndex globalNumbering
    (
        addr.size(),
        Pstream::msgType(),
        mesh.comm(),
        Pstream::parRun()
    );

    const labelListList assemblyCellCells
    (
        assemblyLduMesh.globalCellCells
        (
            assemblyLduMesh,
            globalNumbering
        )
    );

    vectorField cellCentres(addr.size(), Zero);

    forAll(meshes, i)
    {
        const label cellOffset = assemblyLduMesh.cellOffsets()[i];

        forAll(meshes[i].C(), localCellI)
        {
            cellCentres[cellOffset + localCellI] = meshes[i].C()[localCellI];
        }
    }

    labelList assemblyDecomp
    (
        method_().decompose(assemblyCellCells, cellCentres)
    );

    const labelUList& nbr = addr.upperAddr();
    const labelUList& own = addr.lowerAddr();

    bitSet isMappedFace(addr.upperAddr().size());

    forAll(meshes, meshi)
    {
        forAll(meshes[meshi].boundaryMesh(), patchI)
        {
            const polyPatch& pp = meshes[meshi].boundaryMesh()[patchI];

            if (isA<mappedPatchBase>(pp))
            {
                forAll(pp, localFaceI)
                {
                    label allFacei =
                        assemblyLduMesh.faceBoundMap()[meshi][patchI]
                        [
                            localFaceI
                        ];
                    isMappedFace[allFacei] = true;
                }
            }
        }
    }

    while (true)
    {
        label nChanged = 0;

        forAll(own, facei)
        {
            if (isMappedFace[facei])
            {
                label ownProc = assemblyDecomp[own[facei]];
                label neiProc = assemblyDecomp[nbr[facei]];
                if (ownProc < neiProc)
                {
                    assemblyDecomp[nbr[facei]] = ownProc;
                    ++nChanged;
                }
                else if (neiProc < ownProc)
                {
                    assemblyDecomp[own[facei]] = neiProc;
                    ++nChanged;
                }
            }
        }

        if (returnReduce(nChanged, sumOp<label>()) == 0)
        {
            break;
        }
    }

    label offset(0);

    labelList localDecomp(mesh.nCells());

    forAll(meshes, i)
    {
        // Write decompositions for other regions
        if (mesh.name() != meshes[i].name())
        {
            labelIOList otherDecomp
            (
                IOobject
                (
                    dataFile_,
                    meshes[i].facesInstance(),
                    meshes[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                SubList<label>
                (
                    assemblyDecomp,
                    meshes[i].nCells(),
                    offset
                )
            );

            otherDecomp.write();
        }
        else
        {
            localDecomp = SubList<label>(assemblyDecomp, mesh.nCells(), offset);
        }

        offset += meshes[i].nCells();
    }

    // Return decomposition of my region
    return localDecomp;
}


// ************************************************************************* //
