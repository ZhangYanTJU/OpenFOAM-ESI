/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::monitorMappedFilterInletFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    Scalar (e.g. CO2) inlet boundary condition that
     - (area) averages a given scalar field
     (e.g. CO2) over a given "monitor" patch,
     - reduces the average value according to
     a given portion (e.g. 60% CO2 filtering),
     - performs a sharing the resulting value among other
     \c monitorMappedFilterInlet boundary-condition patches, if requested,
     - and applies the resulting value over the inlet patch uniformly
     (e.g. reintroduction of the unfiltered 40% CO2 among various exhausts).

    \f[
        T_{inlet} = T_{monitor} (1 - f_T) \sigma_n
    \f]

    with a constraint for the distribution among patches:
    \f[
        \sum_n \sigma_n = 1
    \f]

    where
    \vartable
      T_{inlet}   | Scalar mass/volume at the inlet patch (area average)
      T_{monitor} | Scalar mass/volume at the monitor patch (area average)
      f_T         | Filtering ratio (i.e. removed)
      \sigma_n    | Share ratio for the patch n
      n           | Patch index among the monitorMappedFilterInlet patches
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type            monitorMappedFilterInlet;
        monitorPatch    <monitorPatchName>;
        filterRatio     0.1;

        // Optional entries (unmodifiable)
        shareRatio      0.5;
        phi             phi;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                            | Type   | Reqd  | Dflt
      type     | Type name: monitorMappedFilterInlet    | word   | yes   | -
      monitorPatch | Name of patch to monitor           | word   | yes   | -
      filterRatio  | Filtering ratio (i.e. removed)     | scalar | no    | 0
      shareRatio   | Share ratio for this patch         | scalar | no    | 1
      phi      | Name of flux field                     | word   | no    | phi
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink

Note
    - Extrema of \c filterRatio and \c shareRatio are \c min=0 and \c max=1.
    - The total \c shareRatio across all \c monitorMappedFilterInlet must be 1.
    - Caveat: No internal check is performed whether the total \c shareRatio
    across all \c monitorMappedFilterInlet equates to 1.

See also
    Foam::fixedValueFvPatchField
    Foam::outletMappedUniformInletHeatAdditionFvPatchField
    Foam::outletMappedUniformInletTemperatureFvPatchField

SourceFiles
    monitorMappedFilterInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef monitorMappedFilterInletFvPatchScalarField_H
#define monitorMappedFilterInletFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class monitorMappedFilterInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class monitorMappedFilterInletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Name of the patch to be monitored
        word monitorPatchName_;

        //- Name of flux field
        word phiName_;

        //- Filtering ratio (i.e. removed)
        scalar filterRatio_;

        //- Share ratio for this patch
        scalar shareRatio_;


public:

    //- Runtime type information
    TypeName("monitorMappedFilterInlet");


    // Constructors

        //- Construct from patch and internal field
        monitorMappedFilterInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        monitorMappedFilterInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- monitorMappedFilterInletFvPatchScalarField
        //- onto a new patch
        monitorMappedFilterInletFvPatchScalarField
        (
            const monitorMappedFilterInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        monitorMappedFilterInletFvPatchScalarField
        (
            const monitorMappedFilterInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new monitorMappedFilterInletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        monitorMappedFilterInletFvPatchScalarField
        (
            const monitorMappedFilterInletFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new monitorMappedFilterInletFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Share ratio for this patch
            const scalar& shareRatio() const
            {
                return shareRatio_;
            }


        // Evaluate

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
