/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) Wikki Ltd
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::pointVolInterpolation

Description

SourceFiles
    pointVolInterpolation.C
    pointVolInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointVolInterpolation_H
#define Foam_pointVolInterpolation_H

#include "primitiveFieldsFwd.H"
#include "primitivePatchInterpolation.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "scalarList.H"
#include "tmp.H"
#include "className.H"
#include "FieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                    Class pointVolInterpolation Declaration
\*---------------------------------------------------------------------------*/

class pointVolInterpolation
{
    // Private Data

        const pointMesh& pointMesh_;

        const fvMesh& fvMesh_;

        //- Interpolation scheme weighting factor array.
        mutable std::unique_ptr<FieldField<Field, scalar>> volWeightsPtr_;

        //- Primitive patch interpolators
        mutable PtrList<primitivePatchInterpolation> patchInterpolators_;


    // Private Member Functions

        //- Return patch interpolators
        const PtrList<primitivePatchInterpolation>& patchInterpolators() const;

        //- Construct point weighting factors
        void makeWeights() const;

        //- Clear addressing
        void clearAddressing() const;

        //- Clear geometry
        void clearGeom() const;


protected:

        const pointMesh& pMesh() const noexcept { return pointMesh_; }

        const fvMesh& vMesh() const noexcept { return fvMesh_; }


public:

    // Declare name of the class and it's debug switch
    ClassName("pointVolInterpolation");


    // Constructors

        //- Construct given pointMesh and fvMesh.
        pointVolInterpolation(const pointMesh&, const fvMesh&);


    //- Destructor
    ~pointVolInterpolation();


    // Member Functions

        // Access

            //- Return reference to weights arrays.
            //  This also constructs the weighting factors if necessary.
            const FieldField<Field, scalar>& volWeights() const;


        // Edit

            //- Update mesh topology using the morph engine
            void updateTopology();

            //- Correct weighting factors for moving mesh.
            bool movePoints();


    // Interpolation functions

        //- Interpolate from pointField to volField
        //  using inverse distance weighting
        template<class Type>
        void interpolate
        (
            const GeometricField<Type, pointPatchField, pointMesh>&,
            GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Interpolate pointField returning volField
        //  using inverse distance weighting
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> interpolate
        (
            const GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate tmp<pointField> returning volField
        //  using inverse distance weighting
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> interpolate
        (
            const tmp<GeometricField<Type, pointPatchField, pointMesh>>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pointVolInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
