/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::pressureInletOutletParSlipVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This velocity inlet/outlet boundary condition for pressure boundary where
    the pressure is specified.  A zero-gradient is applied for outflow (as
    defined by the flux); for inflow, the velocity is obtained from the flux
    with the specified inlet direction.

    A slip condition is applied tangential to the patch.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            pressureInletOutletParSlipVelocity;
        value           uniform 0;
    }
    \endverbatim

Note
    Sign conventions:
    - positive flux (out of domain): apply zero-gradient condition
    - negative flux (into of domain): derive from the flux with specified
      direction

See also
    Foam::mixedFvPatchVectorField
    Foam::pressureDirectedInletOutletVelocityFvPatchVectorField

SourceFiles
    pressureInletOutletParSlipVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureInletOutletParSlipVelocityFvPatchVectorField_H
#define pressureInletOutletParSlipVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class pressureInletOutletParSlipVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pressureInletOutletParSlipVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

        //- Flux field name
        word phiName_;

        //- Density field name
        word rhoName_;


public:

    //- Runtime type information
    TypeName("pressureInletOutletParSlipVelocity");


    // Constructors

        //- Construct from patch and internal field
        pressureInletOutletParSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct by mapping given
        //  pressureInletOutletParSlipVelocityFvPatchVectorField
        //  onto a new patch
        pressureInletOutletParSlipVelocityFvPatchVectorField
        (
            const pressureInletOutletParSlipVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct from patch, internal field and dictionary
        pressureInletOutletParSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct as copy
        pressureInletOutletParSlipVelocityFvPatchVectorField
        (
            const pressureInletOutletParSlipVelocityFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        pressureInletOutletParSlipVelocityFvPatchVectorField
        (
            const pressureInletOutletParSlipVelocityFvPatchVectorField&,
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


    // Member Functions

        //- True: this patch field is altered by assignment
        virtual bool assignable() const { return true; }


        // Access

            //- Return the name of rho
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of rho to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
