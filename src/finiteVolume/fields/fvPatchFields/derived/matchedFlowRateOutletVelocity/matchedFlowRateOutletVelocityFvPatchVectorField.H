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
    Foam::matchedFlowRateOutletVelocityFvPatchVectorField

Group
    grpOutletBoundaryConditions

Description
    Velocity outlet boundary condition which corrects the extrapolated velocity
    to match the flow rate of the specified corresponding inlet patch.

Usage
    \table
        Property     | Description             | Required    | Default value
        inletPatch   | Corresponding inlet patch name | yes  |
        volumetric   | Set volumetric or mass flow-rate | no | false
        rho          | density field name      | no          | rho
    \endtable

    Example of the boundary condition specification for a volumetric flow rate:
    \verbatim
    <patchName>
    {
        type                matchedFlowRateOutletVelocity;
        inletPatch          inlet;
        value               uniform (0 0 0);
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    matchedFlowRateOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef matchedFlowRateOutletVelocityFvPatchVectorField_H
#define matchedFlowRateOutletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class matchedFlowRateOutletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class matchedFlowRateOutletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet patch name from which the corresponding flow rate is obtained
        word inletPatchName_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;

        //- Is volumetric?
        bool volumetric_;


    // Private member functions

        //- Update the patch values given the appropriate density type and value
        template<class RhoType>
        void updateValues
        (
            const label inletPatchID,
            const RhoType& rhoOutlet,
            const RhoType& rhoInlet
        );


public:

   //- Runtime type information
   TypeName("matchedFlowRateOutletVelocity");


   // Constructors

        //- Construct from patch and internal field
        matchedFlowRateOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        matchedFlowRateOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  matchedFlowRateOutletVelocityFvPatchVectorField
        //  onto a new patch
        matchedFlowRateOutletVelocityFvPatchVectorField
        (
            const matchedFlowRateOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        matchedFlowRateOutletVelocityFvPatchVectorField
        (
            const matchedFlowRateOutletVelocityFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        matchedFlowRateOutletVelocityFvPatchVectorField
        (
            const matchedFlowRateOutletVelocityFvPatchVectorField&,
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
