/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

#ifndef Foam_MapConsistentVolFields_H
#define Foam_MapConsistentVolFields_H

#include "GeometricField.H"
#include "meshToMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void evaluateConstraintTypes
(
    GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    fld.boundaryFieldRef().evaluate_if
    (
        [](const auto& pfld) -> bool
        {
            return
            (
                pfld.type() == pfld.patch().patch().type()
             && polyPatch::constraintType(pfld.patch().patch().type())
            );
        },
        UPstream::defaultCommsType
    );
}


template<class Type, class CombineOp>
void MapVolFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields,
    const meshToMesh& interp,
    const CombineOp& cop
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& meshSource = static_cast<const fvMesh&>(interp.srcRegion());
    const fvMesh& meshTarget = static_cast<const fvMesh&>(interp.tgtRegion());

    // Available fields, sorted order
    for
    (
        const IOobject& io :
        (
            selectedFields.empty()
          ? objects.csorted<fieldType>()
          : objects.csorted<fieldType>(selectedFields)
        )
    )
    {
        const fieldType fieldSource(io, meshSource, false);

        IOobject targetIO
        (
            io.name(),
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ
        );

        if (targetIO.typeHeaderOk<fieldType>(true))
        {
            Info<< "    interpolating onto existing field "
                << targetIO.name() << endl;

            fieldType fieldTarget(targetIO, meshTarget, false);

            interp.mapSrcToTgt(fieldSource, cop, fieldTarget);

            evaluateConstraintTypes(fieldTarget);

            fieldTarget.write();
        }
        else
        {
            Info<< "    creating new field "
                << targetIO.name() << endl;

            targetIO.readOpt(IOobject::NO_READ);

            tmp<fieldType> tfieldTarget
            (
                interp.mapSrcToTgt(fieldSource, cop)
            );

            fieldType fieldTarget(targetIO, tfieldTarget);

            evaluateConstraintTypes(fieldTarget);

            fieldTarget.write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
