/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "assemblyFaceAreaPairGAMGAgglomeration.H"
#include "lduMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrixAssembly.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(assemblyFaceAreaPairGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        assemblyFaceAreaPairGAMGAgglomeration,
        lduMatrix
    );

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        assemblyFaceAreaPairGAMGAgglomeration,
        geometry
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::assemblyFaceAreaPairGAMGAgglomeration::
assemblyFaceAreaPairGAMGAgglomeration
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(matrix.mesh(), controlDict)
{
    const fvMatrixAssembly& assembleMatrix =
        static_cast<const fvMatrixAssembly&>(matrix);

    agglomerate
    (
        matrix.mesh(),
        mag
        (
            cmptMultiply
            (
                assembleMatrix.Sf()/sqrt(mag(assembleMatrix.Sf())),
                vector(1, 1.01, 1.02)
            )
        )
    );
}


Foam::assemblyFaceAreaPairGAMGAgglomeration::
assemblyFaceAreaPairGAMGAgglomeration
(
    const lduMatrix& matrix,
    const scalarField& cellVolumes,
    const vectorField& faceAreas,
    const dictionary& controlDict
)
:
    pairGAMGAgglomeration(matrix.mesh(), controlDict)
{
    agglomerate
    (
        matrix.mesh(),
        mag
        (
            cmptMultiply
            (
                faceAreas/sqrt(mag(faceAreas)),
                vector(1, 1.01, 1.02)
            )
        )
    );
}


// ************************************************************************* //
