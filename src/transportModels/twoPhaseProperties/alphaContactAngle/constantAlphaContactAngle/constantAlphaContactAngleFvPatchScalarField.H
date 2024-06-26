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
    Foam::constantAlphaContactAngleFvPatchScalarField

Description
    A constant alphaContactAngle scalar boundary condition.

See also
    Foam::alphaContactAngleTwoPhaseFvPatchScalarField
    Foam::temperatureDependentAlphaContactAngleFvPatchScalarField

SourceFiles
    constantAlphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef constantAlphaContactAngleFvPatchScalarField_H
#define constantAlphaContactAngleFvPatchScalarField_H

#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class constantAlphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class constantAlphaContactAngleFvPatchScalarField
:
    public alphaContactAngleTwoPhaseFvPatchScalarField
{
    // Private data

        //- Equilibrium contact angle
        scalar theta0_;


public:

    //- Runtime type information
    TypeName("constantAlphaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        constantAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        constantAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  constantAlphaContactAngleFvPatchScalarField
        //  onto a new patch
        constantAlphaContactAngleFvPatchScalarField
        (
            const constantAlphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        constantAlphaContactAngleFvPatchScalarField
        (
            const constantAlphaContactAngleFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        constantAlphaContactAngleFvPatchScalarField
        (
            const constantAlphaContactAngleFvPatchScalarField&,
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

        //- Return the equilibrium contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
