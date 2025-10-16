/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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
    Foam::cornerDetectionModels::fluxBased

Description
    Flux-based corner detection model. Marks faces at which the liquid film can
    separate based on the flux.

    The model identifies sharp edges based on the face normals of the two
    faces sharing the edge. Then, if the edge is sharp, the flux direction is
    evaluated to mark the face through which the flux leaves the liquid film.
    If the edge is sharp and the flux leaves through one of the two faces
    sharing the edge, the face is marked as a corner face, where the film can
    separate.

Usage
    Minimal example in boundary-condition files:
    \verbatim
    filmSeparationCoeffs
    {
        // Inherited entries
        ...

        // Optional entries
        cornerDetectionModel  fluxBased;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                     | Type | Reqd | Deflt
      cornerDetectionModel | Corner detector model   | word | no   | fluxBased
    \endtable

    The inherited entries are elaborated in:
      - \link cornerDetectionModel.H \endlink

SourceFiles
    fluxBased.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cornerDetectionModels_fluxBased_H
#define Foam_cornerDetectionModels_fluxBased_H

#include "cornerDetectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cornerDetectionModels
{

/*---------------------------------------------------------------------------*\
                          Class fluxBased Declaration
\*---------------------------------------------------------------------------*/

class fluxBased
:
    public cornerDetectionModel
{
    // Private Data

        //- Flag to deduce if the object is initialised
        bool init_;

        //- Identified corner edges
        bitSet cornerEdges_;


    // Private Member Functions

        //- Return Boolean list of identified corner edges
        bitSet identifyCornerEdges() const;

        //- Return Boolean list of identified corner faces
        bitSet identifyCornerFaces(const bitSet& cornerEdges) const;

        //- Return the list of corner angles for each edge [rad]
        scalarList calcCornerAngles
        (
            const bitSet& faces,
            const bitSet& edges
        ) const;

        // Write edges and faces as OBJ sets for debug purposes
        void writeEdgesAndFaces(const word& prefix = "geomCorners") const;


public:

    //- Runtime type information
    TypeName("fluxBased");


    // Constructors

        //- Construct from components
        fluxBased
        (
            const faMesh& mesh,
            const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
            const dictionary& dict
        );


    // Destructor
    virtual ~fluxBased() = default;


    // Member Functions

    // Evaluation

        //- Detect and store the corner faces
        virtual bool detectCorners();


    // I-O

        //- Read the model dictionary
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cornerDetectionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
