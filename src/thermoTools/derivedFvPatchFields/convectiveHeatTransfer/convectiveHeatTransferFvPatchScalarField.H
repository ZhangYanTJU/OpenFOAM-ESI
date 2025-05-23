/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::compressible::convectiveHeatTransferFvPatchScalarField

Group
    grpCmpBoundaryConditions

Description
    This boundary condition provides a convective heat transfer coefficient
    condition

    if Re > 500000
    \f[
        htc_p = \frac{0.664 Re^{0.5} Pr^{0.333} \kappa_p}{L}
    \f]
    else
    \f[
        htc_p = \frac{0.037 Re^{0.8} Pr^{0.333} \kappa_p}{L}
    \f]

    where

    \vartable
        htc_p   | patch convective heat transfer coefficient
        Re      | Reynolds number
        Pr      | Prandtl number
        \kappa_p | thermal conductivity
        L       | length scale
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        L            | Length scale [m]        | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            convectiveHeatTransfer;
        L               0.1;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    convectiveHeatTransferFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef convectiveHeatTransferFvPatchScalarField_H
#define convectiveHeatTransferFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
           Class convectiveHeatTransferFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class convectiveHeatTransferFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected data

        //- L Length scale [m]
        const scalar L_;


public:

    //- Runtime type information
    TypeName("convectiveHeatTransfer");


    // Constructors

        //- Construct from patch and internal field
        convectiveHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        convectiveHeatTransferFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  convectiveHeatTransferFvPatchScalarField
        //  onto a new patch
        convectiveHeatTransferFvPatchScalarField
        (
            const convectiveHeatTransferFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        convectiveHeatTransferFvPatchScalarField
        (
            const convectiveHeatTransferFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        convectiveHeatTransferFvPatchScalarField
        (
            const convectiveHeatTransferFvPatchScalarField&,
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


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
