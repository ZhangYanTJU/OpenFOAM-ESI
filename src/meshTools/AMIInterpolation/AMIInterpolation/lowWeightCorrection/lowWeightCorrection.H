/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::lowWeightCorrection

Description
    Arbitrary Mesh Interface (AMI) method wrapper that extends low-weight
    faces to use donors of local neighbouring (non-low weight) faces.

    It collects valid nearby faces and uses their stencil instead.

    Example of the boundary specification:
    \verbatim
    <patchName>
    {
        type                cyclicAMI;

        AMIMethod           lowWeightCorrection;

        // Low weight value; which faces to override
        lowWeightCorrection 0.1;
        // How far to walk out from valid faces. 0 means no change.
        nIters              2;
        // Limit amount of valid faces to store.
        nDonors             2;
        // Underlying AMI method to use.
        lowWeightCorrectionCoeffs
        {
            AMIMethod   faceAreaWeightAMI;
        }
    }
    \endverbatim

    This version walks out the local faces and then collects all their
    remote donors. Hence even for e.g. 2 local faces (nDonors) the stencil
    might still include lots of remote faces.

    The number of iterations to walk out might be set to infinite if one
    wants a value on all faces.


SourceFiles
    lowWeightCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef lowWeightCorrection_H
#define lowWeightCorrection_H

#include "AMIInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class lowWeightCorrection Declaration
\*---------------------------------------------------------------------------*/

class lowWeightCorrection
:
    public AMIInterpolation
{
private:

    // Private Data

        //- Threshold weight below which interpolation switches to distance
        //  weighted
        const scalar lowWeightCorrection_;

        //- Number of donors to weigh
        const label nDonors_;

        //- Number of iterations to walk out to search for donors
        const label nIters_;

        //- Uncorrected AMI
        autoPtr<AMIInterpolation> AMIPtr_;


    // Private Member Functions

        //- Find (local) donor faces for all uncovered faces
        void findNearest
        (
            const polyMesh& mesh,
            const primitivePatch& pp,
            const globalIndex& globalFaces,
            const labelList& uncoveredFaces,
            labelListList& nearestCoveredFaces
        ) const;

        //- Determine stencil for single face. Gets list of (local) faces
        //  whose addressing can be used for calculating the stencil
        void calculateStencil
        (
            const label facei,
            const point& sample,
            const labelList& nearestCovered,
            const labelListList& stencilAddressing,

            labelListList& addressing,
            scalarListList& weights,
            scalarField& sumWeights,
            const pointField& faceCentres
        ) const;


        //- No copy assignment
        void operator=(const lowWeightCorrection&) = delete;


public:

    //- Runtime type information
    TypeName("lowWeightCorrection");


    // Constructors

        //- Construct from dictionary
        lowWeightCorrection
        (
            const dictionary& dict,
            const bool reverseTarget = false
        );

        //- Construct from components
        lowWeightCorrection
        (
            const bool requireMatch = true,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1
        );

        //- Construct as copy
        lowWeightCorrection(const lowWeightCorrection& ami);

        //- Construct and return a clone
        virtual autoPtr<AMIInterpolation> clone() const
        {
            return autoPtr<AMIInterpolation>(new lowWeightCorrection(*this));
        }


    //- Destructor
    virtual ~lowWeightCorrection() = default;


    // Member Functions

        //- Update addressing and weights
        virtual bool calculate
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const autoPtr<searchableSurface>& surfPtr = nullptr
        )
        {
            NotImplemented;
            return false;
        }

        //- Update addressing and weights
        virtual bool calculate
        (
            const polyMesh& mesh,
            const label srcPatchi,
            const primitivePatch& srcPatch,
            const label tgtPatchi,
            const primitivePatch& tgtPatch,
            const autoPtr<searchableSurface>& surfPtr = nullptr
        );


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
