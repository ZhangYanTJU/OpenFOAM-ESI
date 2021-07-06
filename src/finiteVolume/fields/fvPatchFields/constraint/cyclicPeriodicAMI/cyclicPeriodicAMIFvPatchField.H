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
    Foam::cyclicPeriodicAMIFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    Version of cyclicAMI. Automatically duplicates itself to obtain full
    coverage.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            cyclicPeriodicAMI;
    }
    \endverbatim

Note
    The outer boundary of the patch pairs must be similar, i.e. if the owner
    patch is transformed to the neighbour patch, the outer perimiter of each
    patch should be identical (or very similar).

See also
    Foam::cyclicAMIFvPatchField
    Foam::AMIInterpolation

SourceFiles
    cyclicPeriodicAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPeriodicAMIFvPatchField_H
#define cyclicPeriodicAMIFvPatchField_H

#include "coupledFvPatchField.H"
#include "cyclicPeriodicAMILduInterfaceField.H"
#include "cyclicPeriodicAMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class cyclicPeriodicAMIFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicPeriodicAMIFvPatchField
:
    virtual public cyclicPeriodicAMILduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const cyclicPeriodicAMIFvPatch& cyclicPeriodicAMIPatch_;


public:

    //- Runtime type information
    TypeName(cyclicPeriodicAMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicPeriodicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicPeriodicAMIFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cyclicPeriodicAMIFvPatchField onto a new
        //  patch
        cyclicPeriodicAMIFvPatchField
        (
            const cyclicPeriodicAMIFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cyclicPeriodicAMIFvPatchField
        (
            const cyclicPeriodicAMIFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new cyclicPeriodicAMIFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cyclicPeriodicAMIFvPatchField
        (
            const cyclicPeriodicAMIFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new cyclicPeriodicAMIFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return local reference cast into the fv patch
            const cyclicPeriodicAMIFvPatch& cyclicPeriodicAMIPatch() const
            {
                return cyclicPeriodicAMIPatch_;
            }


        // Evaluation functions

//            //- Return true if coupled. Note that the underlying patch
//            //  is not coupled() - the points don't align.
//            virtual bool coupled() const;

            //- Return neighbour coupled internal cell data
            virtual tmp<Field<Type>> patchNeighbourField() const;

            //- Return reference to neighbour patchField
            const cyclicPeriodicAMIFvPatchField<Type>& neighbourPatchField()
            const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const bool add,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const;


        // Cyclic AMI coupled interface functions

            //- Does the patch field perform the transformation
            virtual bool doTransform() const
            {
                return
                   !(
                        cyclicPeriodicAMIPatch_.parallel()
                     || pTraits<Type>::rank == 0
                    );
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicPeriodicAMIPatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicPeriodicAMIPatch_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<Type>::rank;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cyclicPeriodicAMIFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
