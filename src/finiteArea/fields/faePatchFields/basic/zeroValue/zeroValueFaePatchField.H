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
    Foam::zeroValueFaePatchField

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
    zeroValueFaePatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zeroValueFaePatchField_H
#define Foam_zeroValueFaePatchField_H

#include "fixedValueFaePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class zeroValueFaePatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class zeroValueFaePatchField
:
    public fixedValueFaePatchField<Type>
{
    //- The parent boundary condition type
    typedef fixedValueFaePatchField<Type> parent_bctype;

public:

    //- Runtime type information
    TypeName("zeroValue");


    // Constructors

        //- Construct from patch and internal field
        zeroValueFaePatchField
        (
            const faPatch& p,
            const DimensionedField<Type, edgeMesh>& iF
        );

        //- Construct from patch, internal field and dictionary
        zeroValueFaePatchField
        (
            const faPatch& p,
            const DimensionedField<Type, edgeMesh>& iF,
            const dictionary& dict
        );

        //- Copy construct onto a new patch (no mapping needed)
        //- onto a new patch
        zeroValueFaePatchField
        (
            const zeroValueFaePatchField<Type>& pfld,
            const faPatch& p,
            const DimensionedField<Type, edgeMesh>& iF,
            const faPatchFieldMapper&
        );

        //- Copy construct with internal field reference
        zeroValueFaePatchField
        (
            const zeroValueFaePatchField<Type>& pfld,
            const DimensionedField<Type, edgeMesh>& iF
        );

        //- Copy construct
        zeroValueFaePatchField(const zeroValueFaePatchField<Type>& pfld)
        :
            zeroValueFaePatchField<Type>(pfld, pfld.internalField())
        {}


        //- Return clone
        virtual tmp<faePatchField<Type>> clone() const
        {
            return faePatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<faePatchField<Type>> clone
        (
            const DimensionedField<Type, edgeMesh>& iF
        ) const
        {
            return faePatchField<Type>::Clone(*this, iF);
        }


    //- Destructor
    virtual ~zeroValueFaePatchField() = default;


    // Member Functions

        //- True: this patch field fixes a value.
        virtual bool fixesValue() const { return true; }

        //- False: this patch field is not altered by assignment.
        virtual bool assignable() const { return false; }

        //- Write (without "value" entry)
        virtual void write(Ostream& os) const;


    // Member Operators

        // Disable assignment operators
        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const faePatchField<Type>&) {}
        virtual void operator+=(const faePatchField<Type>&) {}
        virtual void operator-=(const faePatchField<Type>&) {}
        virtual void operator*=(const faePatchField<scalar>&) {}
        virtual void operator/=(const faePatchField<scalar>&) {}

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
        virtual void operator==(const faePatchField<Type>&) {}
        virtual void operator==(const Field<Type>&) {}
        virtual void operator==(const Type&) {}

        // Prevent automatic comparison rewriting (c++20)
        bool operator!=(const faePatchField<Type>&) const = delete;
        bool operator!=(const Field<Type>&) const = delete;
        bool operator!=(const Type&) const = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zeroValueFaePatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
