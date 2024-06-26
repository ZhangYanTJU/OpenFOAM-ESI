/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
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
    Foam::externalCoupledMixedFvPatchField

Group
    grpGenericBoundaryConditions grpCoupledBoundaryConditions

Description
    Extends the mixed boundary condition with serialisation functions.

    The serialisation functions:

    - writeHeader
    - writeData
    - readData

    It is used for coupling to external applications in combination
    with the externalCoupled functionObject. The default output is one
    line per face, with columns
    \verbatim
        <value> <snGrad> <refValue> <refGrad> <valueFraction>
    \endverbatim

Note
    readData, writeData are not callbacks for regIOobject (since fvPatchField
    not derived from it). They do however do exactly the same - streaming of
    data.

SeeAlso
    mixedFvPatchField
    externalCoupledFunctionObject

SourceFiles
    externalCoupledMixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef externalCoupledMixedFvPatchField_H
#define externalCoupledMixedFvPatchField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class externalCoupledMixedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class externalCoupledMixedFvPatchField
:
    public mixedFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("externalCoupled");


    // Constructors

        //- Construct from patch and internal field
        externalCoupledMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        externalCoupledMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given externalCoupledMixedFvPatchField
        //  onto a new patch
        externalCoupledMixedFvPatchField
        (
            const externalCoupledMixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        externalCoupledMixedFvPatchField
        (
            const externalCoupledMixedFvPatchField&
        );

        //- Construct as copy setting internal field reference
        externalCoupledMixedFvPatchField
        (
            const externalCoupledMixedFvPatchField&,
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


    //- Destructor
    virtual ~externalCoupledMixedFvPatchField() = default;


    // Member functions

        //- Write header
        virtual void writeHeader(Ostream& os) const;

        //- Write data
        virtual void writeData(Ostream& os) const;

        //- Read data
        virtual void readData(Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "externalCoupledMixedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
