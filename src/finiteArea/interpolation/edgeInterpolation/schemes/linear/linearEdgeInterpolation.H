/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::linearEdgeInterpolation

Description
    Central-differencing interpolation scheme class

SourceFiles
    linearEdgeInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef linearEdgeInterpolation_H
#define linearEdgeInterpolation_H

#include "edgeInterpolationScheme.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class linearEdgeInterpolation
:
    virtual public edgeInterpolationScheme<Type>
{
public:

    //- Runtime type information
    TypeName("linear");


    // Generated Methods

        //- No copy construct
        linearEdgeInterpolation(const linearEdgeInterpolation&) = delete;

        //- No copy assignment
        void operator=(const linearEdgeInterpolation&) = delete;


    // Constructors

        //- Construct from mesh
        linearEdgeInterpolation(const faMesh& mesh)
        :
            edgeInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        linearEdgeInterpolation(const faMesh& mesh, Istream&)
        :
            edgeInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        linearEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField&,
            Istream&
        )
        :
            edgeInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const
        {
            return this->mesh().edgeInterpolation::weights();
        }
};


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh>>
linearEdgeInterpolate(const GeometricField<Type, faPatchField, areaMesh>& vf)
{
    return edgeInterpolationScheme<Type>::interpolate
    (
        vf,
        vf.mesh().edgeInterpolation::weights()
    );
}


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh>>
linearEdgeInterpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tinterp =
        linearEdgeInterpolate(tvf());
    tvf.clear();
    return tinterp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
