/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
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
    Foam::JohnsonJacksonParticleSlipFvPatchVectorField

Description
    Partial slip boundary condition for the particulate velocity.

    References:
    \verbatim
        "Multifluid Eulerian modeling of dense gas–solids fluidized bed
         hydrodynamics: Influence of the dissipation parameters"
        Reuge, N.,
        Chemical Engineering Science
        Volume 63, Issue 22, November 2008, pp. 5540-5551
    \endverbatim

    \verbatim
        "Frictional-collisional constitutive relations for granular materials,
         with application to plane shearing"
        Johnson, P.C., and Jackson, R.,
        Journal of Fluid Mechanics
        Volume 176, March 1987, pp. 67-93
    \endverbatim

SourceFiles
    JohnsonJacksonParticleSlipFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef JohnsonJacksonParticleSlipFvPatchVectorField_H
#define JohnsonJacksonParticleSlipFvPatchVectorField_H

#include "partialSlipFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class JohnsonJacksonParticleSlipFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class JohnsonJacksonParticleSlipFvPatchVectorField
:
    public partialSlipFvPatchVectorField
{
    // Private data

        //- Specularity coefficient
        dimensionedScalar specularityCoefficient_;


public:

    //- Runtime type information
    TypeName("JohnsonJacksonParticleSlip");


    // Constructors

        //- Construct from patch and internal field
        JohnsonJacksonParticleSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        JohnsonJacksonParticleSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        JohnsonJacksonParticleSlipFvPatchVectorField
        (
            const JohnsonJacksonParticleSlipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        JohnsonJacksonParticleSlipFvPatchVectorField
        (
            const JohnsonJacksonParticleSlipFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        JohnsonJacksonParticleSlipFvPatchVectorField
        (
            const JohnsonJacksonParticleSlipFvPatchVectorField&,
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

        //- Update the coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
