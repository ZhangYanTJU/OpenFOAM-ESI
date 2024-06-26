/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::interstitialInletVelocityFvPatchVectorField

Description
    Inlet velocity in which the actual interstitial velocity is calculated
    by dividing the specified inletVelocity field with the local phase-fraction.

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type              interstitialInletVelocity;
        inletVelocity     uniform (0 0.2 0);// Non-interstitial inlet velocity
        alpha             alpha.particles;  // Name of the phase-fraction field
        value             uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    interstitialInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef interstitialInletVelocityFvPatchVectorField_H
#define interstitialInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class interstitialInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class interstitialInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet velocity
        vectorField inletVelocity_;

        //- Name of the flux transporting the field
        word alphaName_;


public:

   //- Runtime type information
   TypeName("interstitialInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        interstitialInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interstitialInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  interstitialInletVelocityFvPatchVectorField
        //  onto a new patch
        interstitialInletVelocityFvPatchVectorField
        (
            const interstitialInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        interstitialInletVelocityFvPatchVectorField
        (
            const interstitialInletVelocityFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        interstitialInletVelocityFvPatchVectorField
        (
            const interstitialInletVelocityFvPatchVectorField&,
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

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


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
