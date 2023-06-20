/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2025 OpenCFD Ltd.
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
    Foam::zeroValueFvsPatchField

Description
    Specifies a zero fixed value boundary condition.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type    zeroValue;
    }
    \endverbatim

SourceFiles
    zeroValueFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zeroValueFvsPatchField_H
#define Foam_zeroValueFvsPatchField_H

#include "fixedValueFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class zeroValueFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zeroValueFvsPatchField
:
    public fixedValueFvsPatchField<Type>
{
    //- The parent boundary condition type
    typedef fixedValueFvsPatchField<Type> parent_bctype;

public:

    //- Runtime type information
    TypeName("zeroValue");


    // Constructors

        //- Construct from patch and internal field
        zeroValueFvsPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, surfaceMesh>& iF
        );

        //- Construct from patch, internal field and dictionary
        zeroValueFvsPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, surfaceMesh>& iF,
            const dictionary& dict
        );

        //- Copy construct onto a new patch (no mapping needed)
        zeroValueFvsPatchField
        (
            const zeroValueFvsPatchField<Type>& pfld,
            const fvPatch& p,
            const DimensionedField<Type, surfaceMesh>& iF,
            const fvPatchFieldMapper&
        );

        //- Copy construct with internal field reference
        zeroValueFvsPatchField
        (
            const zeroValueFvsPatchField<Type>& pfld,
            const DimensionedField<Type, surfaceMesh>& iF
        );

        //- Copy construct
        zeroValueFvsPatchField(const zeroValueFvsPatchField<Type>& pfld)
        :
            zeroValueFvsPatchField<Type>(pfld, pfld.internalField())
        {}


        //- Return clone
        virtual tmp<fvsPatchField<Type>> clone() const
        {
            return fvsPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return fvsPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        //- True: this patch field fixes a value.
        virtual bool fixesValue() const { return true; }

        //- False: this patch field is not altered by assignment.
        virtual bool assignable() const { return false; }

        //- Write (without "value" entry)
        virtual void write(Ostream& os) const;


    // Evaluation Functions

        //- Return the matrix diagonal coefficients corresponding to the
        //- evaluation of the value of this patchField with given weights
        virtual tmp<Field<Type>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //- evaluation of the value of this patchField with given weights
        virtual tmp<Field<Type>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the
        //- evaluation of the gradient of this patchField
        virtual tmp<Field<Type>> gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //- evaluation of the gradient of this patchField
        virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;


    // Member Operators

        // Disable assignment operators
        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvsPatchField<Type>&) {}
        virtual void operator+=(const fvsPatchField<Type>&) {}
        virtual void operator-=(const fvsPatchField<Type>&) {}
        virtual void operator*=(const fvsPatchField<scalar>&) {}
        virtual void operator/=(const fvsPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}


        // Disable forced assignment operators
        virtual void operator==(const fvsPatchField<Type>&) {}
        virtual void operator==(const Field<Type>&) {}
        virtual void operator==(const Type&) {}

        // Prevent automatic comparison rewriting (c++20)
        bool operator!=(const fvPatchField<Type>&) const = delete;
        bool operator!=(const Field<Type>&) const = delete;
        bool operator!=(const Type&) const = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zeroValueFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
