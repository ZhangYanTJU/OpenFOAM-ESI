/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 PCOpt/NTUA
    Copyright (C) 2022      FOSS GP
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

Class
    shortGeometricField

Description
    Stores only the current state of a GeometricField (i.e. not oldTimes).
    Additionally, avoids storing patchFields that can be retrieved through
    correctBoundaryConditions (e.g. uniform fields, coupled, symmetry, etc)

SourceFiles
    shortGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef shortGeometricField_H
#define shortGeometricField_H

#include "compressedGeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class shortGeometricField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class shortGeometricField
:
    public compressedGeometricField<Type, PatchField, GeoMesh>
{
protected:

    // Protected data

        //- Storage of the internalField deceomposed into multiple
        //- scalarFields
        List<scalarField> internalBufferPtr_;

        //- Storage of the boundaryField deceomposed into multiple
        //- scalarFields per patch
        List<List<scalarField>> boundaryBufferPtr_;

        //- Holds the uniform value for patches with such boundary conditions.
        //  Allocated only on these patches
        PtrList<Type> uniformBoundaryValuePtr_;


    // Protected Member Functions

        //- Fill in uniformBoundaryValuePtr_ with the uniform boundary field
        //- values
        void storeUniformBoundaryValues() override;

        //- Restore the uniform boundary field values based on what is stored
        //- in uniformBoundaryValuePtr_
        void restoreUniformBoundaryValues
        (
            GeometricField<Type, PatchField, GeoMesh>& field
        ) override;

        //- Initialize necessary parameters
        void initialize();

        //- Set sizes for internalBufferPtr_ and boundaryBufferPtr_
        void setBuffers();


public:

    //- Runtime type information
    TypeName("shortGeometricField");


    // Constructors
    shortGeometricField
    (
        GeometricField<Type, PatchField, GeoMesh>& field,
        storageParameters& storageParams,
        const label iPtr,
        const label& k
    );


    // Destructors
    virtual ~shortGeometricField() = default;


    // Member Functions

        //- Store field
        void compress() override;

        //- Retrieve field
        void decompress
        (
            GeometricField<Type, PatchField, GeoMesh>& field
        ) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "shortGeometricField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
