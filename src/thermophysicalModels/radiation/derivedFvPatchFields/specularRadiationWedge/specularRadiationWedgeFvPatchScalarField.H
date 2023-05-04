/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Dominik Christ, Wikki Ltd.
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::specularRadiationWedgeFvPatchScalarField

Description
    Foam::specularRadiationWedgeFvPatchScalarField

SourceFiles
    specularRadiationWedgeFvPatchScalarField.C

Author
    Dominik Christ, Wikki Ltd. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef specularRadiationWedgeFvPatchScalarField_H
#define specularRadiationWedgeFvPatchScalarField_H

#include "wedgeFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
            Class specularRadiationWedgeFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class specularRadiationWedgeFvPatchScalarField
:
    public wedgeFvPatchField<scalar>
{

    //- For each face, ray index of ray in direction of reflection
    mutable labelListList *receivedRayIDPtr_;

    //- Calculate IDs of rays in direction of reflection
    void calcReceivedRayIDs() const;


public:

    //- Runtime type information
    TypeName("specularRadiationWedge");


    // Constructors

        //- Construct from patch and internal field
        specularRadiationWedgeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        specularRadiationWedgeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given specularRadiationWedgeFvPatchScalarField onto a new patch
        specularRadiationWedgeFvPatchScalarField
        (
            const specularRadiationWedgeFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        specularRadiationWedgeFvPatchScalarField
        (
            const specularRadiationWedgeFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new specularRadiationWedgeFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        specularRadiationWedgeFvPatchScalarField
        (
            const specularRadiationWedgeFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new specularRadiationWedgeFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return ray IDs of reflected rays
            const labelListList& receivedRayIDs() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchScalarField& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
