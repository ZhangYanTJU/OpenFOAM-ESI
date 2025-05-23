/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022-2025 OpenCFD Ltd.
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

#include "fvPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeometricField, class AnyType>
const typename GeometricField::Patch&
Foam::fvPatch::lookupPatchField(const word& fldName) const
{
    return
        boundaryMesh().mesh().thisDb().template
            lookupObject<GeometricField>(fldName)
            .boundaryField()[this->index()];
}


template<class GeometricField>
const typename GeometricField::Patch*
Foam::fvPatch::cfindPatchField(const word& fldName) const
{
    const auto* fldptr =
        boundaryMesh().mesh().thisDb().template
        cfindObject<GeometricField>(fldName);

    if (fldptr)
    {
        return &(fldptr->boundaryField()[this->index()]);
    }
    else
    {
        return nullptr;
    }
}


// ************************************************************************* //
