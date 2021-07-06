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

Class
    Foam::cyclicPeriodicAMIFvPatch

Description
    Periodic cyclic patch for Arbitrary Mesh Interface (AMI)

SourceFiles
    cyclicPeriodicAMIFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPeriodicAMIFvPatch_H
#define cyclicPeriodicAMIFvPatch_H

#include "coupledFvPatch.H"
#include "cyclicPeriodicAMILduInterface.H"
#include "cyclicPeriodicAMIPolyPatch.H"
#include "fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicPeriodicAMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicPeriodicAMIFvPatch
:
    public coupledFvPatch,
    public cyclicPeriodicAMILduInterface
{
    // Private data

        const cyclicPeriodicAMIPolyPatch& cyclicPeriodicAMIPolyPatch_;


protected:

    // Protected Member functions

        //- Note: straight copy of cyclicAMIFvPatch

        //- Make patch weighting factors
        void makeWeights(scalarField&) const;

        //- Correct patch deltaCoeffs
        virtual void makeDeltaCoeffs(scalarField&) const;

        //- Correct patch non-ortho deltaCoeffs
        virtual void makeNonOrthoDeltaCoeffs(scalarField&) const;

        //- Correct patch non-ortho correction vectors
        virtual void makeNonOrthoCorrVectors(vectorField&) const;

        //- Correct patches after moving points
        virtual void movePoints();


public:

    //- Runtime type information
    TypeName(cyclicPeriodicAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        cyclicPeriodicAMIFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            coupledFvPatch(patch, bm),
            cyclicPeriodicAMILduInterface(),
            cyclicPeriodicAMIPolyPatch_
            (
                refCast<const cyclicPeriodicAMIPolyPatch>(patch)
            )
        {}


    // Member functions

        // Access

            //- Return local reference cast into the cyclic patch
            const cyclicPeriodicAMIPolyPatch& cyclicPeriodicAMIPatch() const
            {
                return cyclicPeriodicAMIPolyPatch_;
            }

            //- Return neighbour
            virtual label neighbPatchID() const
            {
                return cyclicPeriodicAMIPolyPatch_.neighbPatchID();
            }

            virtual bool owner() const
            {
                return cyclicPeriodicAMIPolyPatch_.owner();
            }

            //- Return processor number
            virtual const cyclicPeriodicAMIFvPatch& neighbPatch() const
            {
                return refCast<const cyclicPeriodicAMIFvPatch>
                (
                    this->boundaryMesh()
                    [
                        cyclicPeriodicAMIPolyPatch_.neighbPatchID()
                    ]
                );
            }

            //- Return a reference to the AMI interpolator
            virtual const AMIPatchToPatchInterpolation& AMI() const
            {
                return cyclicPeriodicAMIPolyPatch_.AMI();
            }

            //- Return true if applying the low weight correction
            virtual bool applyLowWeightCorrection() const
            {
                return cyclicPeriodicAMIPolyPatch_.applyLowWeightCorrection();
            }

            //- Are the cyclic planes parallel
            virtual bool parallel() const
            {
                return cyclicPeriodicAMIPolyPatch_.parallel();
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicPeriodicAMIPolyPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicPeriodicAMIPolyPatch_.reverseT();
            }

            const cyclicPeriodicAMIFvPatch& neighbFvPatch() const
            {
                return refCast<const cyclicPeriodicAMIFvPatch>
                (
                    this->boundaryMesh()
                    [
                        cyclicPeriodicAMIPolyPatch_.neighbPatchID()
                    ]
                );
            }

            //- Return true if this patch is coupled. This is equivalent
            //- to the coupledPolyPatch::coupled() if parallel running or
            //- both sides present, false otherwise
            virtual bool coupled() const;

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& fld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return cyclicPeriodicAMIPolyPatch_.interpolate
                (
                    fld,
                    defaultValues
                );
            }

            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const tmp<Field<Type>>& tFld,
                const UList<Type>& defaultValues = UList<Type>()
            ) const
            {
                return cyclicPeriodicAMIPolyPatch_.interpolate
                (
                    tFld,
                    defaultValues
                );
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //- the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
