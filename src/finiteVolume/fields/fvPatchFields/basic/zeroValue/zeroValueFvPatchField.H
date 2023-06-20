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
    Foam::zeroValueFvPatchField

Group
    grpGenericBoundaryConditions

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
    zeroValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zeroValueFvPatchField_H
#define Foam_zeroValueFvPatchField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class zeroValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zeroValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    //- The parent boundary condition type
    typedef fixedValueFvPatchField<Type> parent_bctype;

public:

    //- Runtime type information
    TypeName("zeroValue");


    // Constructors

        //- Construct from patch and internal field
        zeroValueFvPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF
        );

        //- Construct from patch, internal field and dictionary
        zeroValueFvPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF,
            const dictionary& dict
        );

        //- Copy construct onto a new patch (no mapping needed)
        zeroValueFvPatchField
        (
            const zeroValueFvPatchField<Type>& pfld,
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF,
            const fvPatchFieldMapper&
        );

        //- Copy construct with internal field reference
        zeroValueFvPatchField
        (
            const zeroValueFvPatchField<Type>& pfld,
            const DimensionedField<Type, volMesh>& iF
        );

        //- Copy construct
        zeroValueFvPatchField(const zeroValueFvPatchField<Type>& pfld)
        :
            zeroValueFvPatchField<Type>(pfld, pfld.internalField())
        {}


        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new zeroValueFvPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new zeroValueFvPatchField<Type>(*this, iF)
            );
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

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

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
        virtual void operator==(const fvPatchField<Type>&) {}
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
    #include "zeroValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
