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

#include "cyclicAMIPointPatchFields.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
//#include "cyclicPeriodicAMIPointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//makePointPatchFields(cyclicPeriodicAMI);


// Redirect cyclicPeriodicAMI to cyclicAMI for now
#define addNamedToPointPatchFieldRunTimeSelection\
(PatchTypeField,typePatchTypeField,lookup)                                     \
    addNamedToRunTimeSelectionTable                                            \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        pointPatch,                                                            \
        lookup                                                                 \
    );                                                                         \
    addNamedToRunTimeSelectionTable                                            \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper,                                                           \
        lookup                                                                 \
    );                                                                         \
    addNamedToRunTimeSelectionTable                                            \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary,                                                            \
        lookup                                                                 \
    );

#define makeNamedPointPatchFields(type,lookupType)                             \
    addNamedToPointPatchFieldRunTimeSelection                                  \
    (                                                                          \
        pointPatchScalarField,                                                 \
        type##PointPatchScalarField,                                           \
        lookupType                                                             \
    );                                                                         \
    addNamedToPointPatchFieldRunTimeSelection                                  \
    (                                                                          \
        pointPatchVectorField,                                                 \
        type##PointPatchVectorField,                                           \
        lookupType                                                             \
    );                                                                         \
    addNamedToPointPatchFieldRunTimeSelection                                  \
    (                                                                          \
        pointPatchSphericalTensorField,                                        \
        type##PointPatchSphericalTensorField,                                  \
        lookupType                                                             \
    );                                                                         \
    addNamedToPointPatchFieldRunTimeSelection                                  \
    (                                                                          \
        pointPatchSymmTensorField,                                             \
        type##PointPatchSymmTensorField,                                       \
        lookupType                                                             \
    );                                                                         \
    addNamedToPointPatchFieldRunTimeSelection                                  \
    (                                                                          \
        pointPatchTensorField,                                                 \
        type##PointPatchTensorField,                                           \
        lookupType                                                             \
);

makeNamedPointPatchFields(cyclicAMI, cyclicPeriodicAMI);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
