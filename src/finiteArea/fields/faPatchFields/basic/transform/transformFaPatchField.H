/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::transformFaPatchField

Description

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    transformFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef transformFaPatchField_H
#define transformFaPatchField_H

#include "faPatchField.H"
#include "areaFaMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class transformFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class transformFaPatchField
:
    public faPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("transform");


    // Constructors

        //- Construct from patch and internal field
        transformFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary.
        //  The "value" entry: NO_READ
        transformFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given transform patch field
        //- onto a new patch
        transformFaPatchField
        (
            const transformFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        transformFaPatchField
        (
            const transformFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return clone
        virtual tmp<faPatchField<Type>> clone() const = 0;

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>&
        ) const = 0;


    // Member functions

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const = 0;

            //- Return face-gradient transform diagonal
            virtual tmp<Field<Type>> snGradTransformDiag() const = 0;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const;


    // Member operators

        virtual void operator=(const faPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "transformFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
