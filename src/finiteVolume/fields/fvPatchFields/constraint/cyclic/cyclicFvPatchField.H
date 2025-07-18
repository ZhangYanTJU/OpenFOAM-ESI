/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::cyclicFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition enforces a cyclic condition between a pair of
    boundaries.

    The "value" entry is NO_READ.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            cyclic;
    }
    \endverbatim

Note
    The patches must be topologically similar, i.e. if the owner patch is
    transformed to the neighbour patch, the patches should be identical (or
    very similar).

SourceFiles
    cyclicFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cyclicFvPatchField_H
#define Foam_cyclicFvPatchField_H

#include "coupledFvPatchField.H"
#include "cyclicLduInterfaceField.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicFvPatchField
:
    virtual public cyclicLduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private Data

        //- Local reference cast into the cyclic patch
        const cyclicFvPatch& cyclicPatch_;


public:


    //- Runtime type information
    TypeName(cyclicFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&,
            const bool needValue = true  //!< evaluate
        );

        //- Construct by mapping given cyclicFvPatchField onto a new patch
        cyclicFvPatchField
        (
            const cyclicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cyclicFvPatchField
        (
            const cyclicFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        cyclicFvPatchField
        (
            const cyclicFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        // Access

            //- Return local reference cast into the cyclic patch
            const cyclicFvPatch& cyclicPatch() const
            {
                return cyclicPatch_;
            }


        // Evaluation functions

            //- Return neighbour coupled internal cell data
            virtual tmp<Field<Type>> patchNeighbourField() const;

            //- Retrieve neighbour coupled internal cell data
            virtual void patchNeighbourField(UList<Type>& pnf) const;

            //- Return reference to neighbour patchField
            const cyclicFvPatchField<Type>& neighbourPatchField() const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const Field<Type>& psiInternal,
                const scalarField& coeffs,
                const Pstream::commsTypes commsType
            ) const;


            //- Manipulate matrix
            virtual void manipulateMatrix
            (
                fvMatrix<Type>& m,
                const label iMatrix,
                const direction cmp
            );



        // Cyclic coupled interface functions

            //- Does the patch field perform the transformation
            virtual bool doTransform() const
            {
                return (pTraits<Type>::rank && !cyclicPatch_.parallel());
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicPatch_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cyclicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
