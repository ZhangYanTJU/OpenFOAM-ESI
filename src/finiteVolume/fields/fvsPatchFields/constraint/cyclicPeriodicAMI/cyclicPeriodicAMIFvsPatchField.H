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
    Foam::cyclicPeriodicAMIFvsPatchField

SourceFiles
    cyclicPeriodicAMIFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicPeriodicAMIFvsPatchField_H
#define cyclicPeriodicAMIFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "cyclicPeriodicAMIFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class cyclicPeriodicAMIFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cyclicPeriodicAMIFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private data

        //- Local reference cast into the cyclic patch
        const cyclicPeriodicAMIFvPatch& cyclicPeriodicAMIPatch_;


public:

    //- Runtime type information
    TypeName(cyclicPeriodicAMIFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        cyclicPeriodicAMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cyclicPeriodicAMIFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cyclicPeriodicAMIFvsPatchField onto
        //  a new patch
        cyclicPeriodicAMIFvsPatchField
        (
            const cyclicPeriodicAMIFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cyclicPeriodicAMIFvsPatchField
        (
            const cyclicPeriodicAMIFvsPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvsPatchField<Type>> clone() const
        {
            return tmp<fvsPatchField<Type>>
            (
                new cyclicPeriodicAMIFvsPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cyclicPeriodicAMIFvsPatchField
        (
            const cyclicPeriodicAMIFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvsPatchField<Type>> clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type>>
            (
                new cyclicPeriodicAMIFvsPatchField<Type>(*this, iF)
            );
        }

    // Member functions

        // Access

            //- Return true if running parallel
            virtual bool coupled() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cyclicPeriodicAMIFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
