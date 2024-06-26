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
    Foam::fluxCorrectedVelocityFvPatchVectorField

Group
    grpOutletBoundaryConditions

Description
    This boundary condition provides a velocity outlet boundary condition for
    patches where the pressure is specified.  The outflow velocity is obtained
    by "zeroGradient" and then corrected from the flux:

        \f[
            U_p = U_c - n (n \cdot U_c) + \frac{n \phi_p}{|S_f|}
        \f]

    where

    \vartable
        U_p | velocity at the patch [m/s]
        U_c | velocity in cells adjacent to the patch [m/s]
        n   | patch normal vectors
        \phi_p | flux at the patch [m3/s or kg/s]
        S_f | patch face area vectors [m2]
    \endvartable

    where

    \table
        Property     | Description             | Required    | Default value
        phi          | name of flux field      | no          | phi
        rho          | name of density field   | no          | rho
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fluxCorrectedVelocity;
        phi             phi;
        rho             rho;
    }
    \endverbatim

Note
    If reverse flow is possible or expected use the
    pressureInletOutletVelocity condition instead.

See also
    Foam::zeroGradientFvPatchField
    Foam::pressureInletOutletVelocityFvPatchVectorField

SourceFiles
    fluxCorrectedVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fluxCorrectedVelocityFvPatchVectorField_H
#define fluxCorrectedVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class fluxCorrectedVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fluxCorrectedVelocityFvPatchVectorField
:
    public zeroGradientFvPatchVectorField
{
    // Private data

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;


public:

    //- Runtime type information
    TypeName("fluxCorrectedVelocity");


    // Constructors

        //- Construct from patch and internal field
        fluxCorrectedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fluxCorrectedVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fluxCorrectedVelocityFvPatchVectorField
        //  onto a new patch
        fluxCorrectedVelocityFvPatchVectorField
        (
            const fluxCorrectedVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        fluxCorrectedVelocityFvPatchVectorField
        (
            const fluxCorrectedVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    // Member functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::buffered
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
