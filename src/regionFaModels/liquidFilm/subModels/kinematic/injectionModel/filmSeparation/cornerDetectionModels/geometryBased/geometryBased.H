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
    Foam::cornerDetectionModels::geometryBased

Description
    Geometry-based corner detection model. Marks faces at which the liquid film
    can separate based on the geometry.

    The model identifies sharp edges based on the face normals of the two
    faces sharing the edge. Then, if the edge is sharp, the curvature is
    evaluated to mark the face through which the flux leaves the liquid film.
    If the edge is sharp and the curvature is of the specified type, the face
    is marked as a corner face, where the film can separate.

Usage
    Minimal example in boundary-condition files:
    \verbatim
    filmSeparationCoeffs
    {
        // Inherited entries
        ...

        // Optional entries
        cornerDetectionModel  geometryBased;
        cornerCurveType       <word>;
        cornerType            <word>;
        angleSharp            <scalar>;   // [deg]
        angleRoundMin         <scalar>;   // [deg]
        angleRoundMax         <scalar>;   // [deg]
        maxRoundRadius        <scalar>;   // [m]
        minEdgeLength         <scalar>;   // [m]
        nSmooth               <label>;    // [no. of passes]
        sharpBoundaryEdges    <bool>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                   | Type | Reqd | Deflt
      cornerDetectionModel | Corner detector model | word | no   | fluxBased
      cornerCurveType | Corner-curvature type      | word | no   | any
      cornerType      | Corner type                | word | no   | sharpOrRound
      angleSharp      | Sharp-angle limit [deg]    | scalar | no | 45
      angleRoundMin   | Minimum round-angle limit [deg] | scalar | no | 5
      angleRoundMax   | Maximum round-angle limit [deg] | scalar | no | 45
      maxRoundRadius  | Maximum round-radius limit [m] | scalar | no | 2e-3
      minEdgeLength   | Minimum edge length [m]     | scalar | no | 0
      nSmooth         | No. of smoothing passes     | label | no | 0
      sharpBoundaryEdges | Treat boundary edges as sharp | bool | no | false
    \endtable

    Options for the \c cornerCurve entry:
    \verbatim
        any          | Convex or concave corners
        convex       | Convex corners
        concave      | Concave corners
    \endverbatim

    Options for the \c cornerType entry:
    \verbatim
        sharp        | Sharp corners
        round        | Round corners
        sharpOrRound | Sharp or round corners
    \endverbatim

    The inherited entries are elaborated in:
      - \link cornerDetectionModel.H \endlink

SourceFiles
    geometryBased.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cornerDetectionModels_geometryBased_H
#define Foam_cornerDetectionModels_geometryBased_H

#include "cornerDetectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cornerDetectionModels
{

/*---------------------------------------------------------------------------*\
                          Class geometryBased Declaration
\*---------------------------------------------------------------------------*/

class geometryBased
:
    public cornerDetectionModel
{
    // Private Enumerations

        //- Options for the corner-curvature type
        enum cornerCurveType : char
        {
            ANY = 0,    //!< "Convex or concave corners"
            CONVEX,     //!< "Convex corners"
            CONCAVE     //!< "Concave corners"
        };

        //- Names for cornerCurveType
        static const Enum<cornerCurveType> cornerCurveTypeNames;


        //- Options for the corner type
        enum cornerType : char
        {
            ALL = 0,    //!< "Sharp or round corners"
            SHARP,      //!< "Sharp corners"
            ROUND       //!< "Round corners"
        };

        //- Names for cornerType
        static const Enum<cornerType> cornerTypeNames;


    // Private Data

        //- Corner-curvature type
        enum cornerCurveType cornerCurveType_;

        //- Corner type
        enum cornerType cornerType_;

        //- Flag to deduce if the object is initialised
        bool init_;

        //- Bitset of edges identified as sharp
        bitSet sharpEdges_;

        //- Bitset of edges identified as round
        bitSet roundEdges_;

        //- Bitset of edges identified as a combination of sharp and round
        bitSet cornerEdges_;

        //- Bitset of faces identified as sharp
        bitSet sharpFaces_;

        //- Bitset of faces identified as round
        bitSet roundFaces_;

        //- Bitset of faces identified as a combination of sharp and round
        bitSet cornerFaces_;

        //- Sharp-angle limit
        scalar angleSharpDeg_;

        //- Minimum round-angle limit
        scalar angleRoundMinDeg_;

        //- Maximum round-angle limit
        scalar angleRoundMaxDeg_;

        //- Maximum round-radius limit
        scalar maxRoundRadius_;

        //- Minimum edge length; ignore edges shorter than this
        scalar minEdgeLength_;

        //- Number of smoothing passes on the binary edge mask
        label nSmooth_;

        //- Flag to treat one-face boundary edges as sharp
        bool sharpBoundaryEdges_;


    // Private Member Functions

        //- Return the signed bending sense, sign(+1/-1) of curvature, across
        //- an edge with respect to the edge tangent
        scalar curvatureSign
        (
            const vector& t,
            const vector& n0,
            const vector& n1
        ) const;

        //- Classify edges into sharp/round according to dihedral angle and
        //- inferred radius
        void classifyEdges
        (
            bitSet& sharpEdges,
            bitSet& roundEdges
        ) const;

        //- Convert an edge mask to a face mask (face is set if any of its
        //- edges are set)
        void edgesToFaces
        (
            const bitSet& edgeMask,
            bitSet& faceMask
        ) const;

        //- Return the list of corner angles [rad] for each edge
        scalarList calcCornerAngles
        (
            const bitSet& faces,
            const bitSet& edges
        ) const;

        // Write edges and faces as OBJ sets for debug purposes
        void writeEdgesAndFaces() const;


public:

    //- Runtime type information
    TypeName("geometryBased");


    // Constructors

        //- Construct from components
        geometryBased
        (
            const faMesh& mesh,
            const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
            const dictionary& dict
        );


    // Destructor
    virtual ~geometryBased() = default;


    // Member Functions

    // Evaluation

        //- Detect and store the corner faces
        virtual bool detectCorners();


    // I-O

        //- Read the model dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cornerDetectionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
