/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::interfaceCompressionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    Applies interface-compression to the phase-fraction distribution at the
    patch by setting the phase-fraction to 0 if it is below 0.5, otherwise
    to 1.

    This approach is useful to avoid unphysical "bleed" of the lighter phase
    along the surface in regions of high shear adjacent to the surface which
    is at a shallow angle to the interface.

See also
    Foam::fixedValueFvPatchScalarField

SourceFiles
    interfaceCompressionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCompressionFvPatchScalarField_H
#define interfaceCompressionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class interfaceCompressionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class interfaceCompressionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("interfaceCompression");


    // Constructors

        //- Construct from patch and internal field
        interfaceCompressionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interfaceCompressionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  interfaceCompressionFvPatchScalarField onto a new patch
        interfaceCompressionFvPatchScalarField
        (
            const interfaceCompressionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        interfaceCompressionFvPatchScalarField
        (
            const interfaceCompressionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        interfaceCompressionFvPatchScalarField
        (
            const interfaceCompressionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
