/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::highAspectRatioFvGeometryScheme

Description
    Default geometry calculation scheme. Slight stabilisation for bad meshes.

SourceFiles
    highAspectRatioFvGeometryScheme.C

\*---------------------------------------------------------------------------*/

#ifndef highAspectRatioFvGeometryScheme_H
#define highAspectRatioFvGeometryScheme_H

#include "basicFvGeometryScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class primitiveMesh;

/*---------------------------------------------------------------------------*\
               Class highAspectRatioFvGeometryScheme Declaration
\*---------------------------------------------------------------------------*/

class highAspectRatioFvGeometryScheme
:
    public basicFvGeometryScheme
{
    // Private Member Functions

        void cellClosedness
        (
            const vectorField& areas,
            const scalarField& vols,
            const tensorField& cellCoords,
            scalarField& aratio
        ) const;

        void cellDirections
        (
            tensorField& T,
            vectorField& lambda
        ) const;


        //- No copy construct
        highAspectRatioFvGeometryScheme
        (
            const highAspectRatioFvGeometryScheme&
        ) = delete;

        //- No copy assignment
        void operator=(const highAspectRatioFvGeometryScheme&) = delete;


public:

    //- Runtime type information
    TypeName("highAspectRatio");


    // Constructors

        //- Construct from mesh
        highAspectRatioFvGeometryScheme
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~highAspectRatioFvGeometryScheme() = default;


    // Member Functions

        //- Do what is necessary if the mesh has moved
        virtual void movePoints();

        //- Helper : calculate (weighted) average face and cell centres
        static void makeAverageCentres
        (
            const primitiveMesh& mesh,
            const pointField& points,
            const pointField& faceAreas,
            pointField& faceCentres,
            pointField& cellCentres
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
