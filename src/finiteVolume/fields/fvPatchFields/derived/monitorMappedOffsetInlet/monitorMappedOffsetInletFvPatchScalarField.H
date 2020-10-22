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
    Foam::monitorMappedOffsetInletFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    Scalar (e.g. temperature) inlet boundary condition that
     - averages the given scalar field over a given "monitor" patch,
     - adds/subtracts a given scalar value unless a given upper/lower
     limit of the scalar is reached,
     - and applies the offset-scalar value over the inlet patch uniformly.

    \f[
        T_{inlet} = T_{monitor} + T_{input}
    \f]

    with the maximal constraint:
    \f[
        T_{inlet} \leq T_{max}
    \f]

    with the minimal constaint:
    \f[
        T_{inlet} \geq T_{min}
    \f]

    where
    \vartable
      T_{inlet}   | Scalar at the inlet patch
      T_{monitor} | Scalar at the monitor patch
      T_{input}   | User-defined scalar to add/subtract
      T_{max}     | Upper-limit scalar (e.g. temperature to cease heating)
      T_{min}     | Lower-limit scalar (e.g. temperature to cease cooling)
    \endvartable

    Required fields:
    \verbatim
      T   | A scalar field (e.g. temperature)
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type            monitorMappedOffsetInlet;
        monitorPatch    <monitorPatchName>;
        offset          10;

        // Optional entries (unmodifiable)
        phi             phi;
        max             1000;
        min             0;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                           | Type   | Reqd  | Dflt
      type     | Type name: monitorMappedOffsetInlet   | word   | yes   | -
      monitorPatch | Name of monitor patch             | word   | yes   | -
      offset   | Scalar value to add or subtract       | scalar | yes   | -
      max      | Upper-limit scalar                    | scalar | no    | GREAT
      min      | Lower-limit scalar                    | scalar | no    | SMALL
      phi      | Name of flux field                    | word   | no    | phi
    \endtable

    The inherited entries are elaborated in:
     - \link fixedValueFvPatchFields.H \endlink

See also
    Foam::fixedValueFvPatchField
    Foam::outletMappedUniformInletHeatAdditionFvPatchField

SourceFiles
    monitorMappedOffsetInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef monitorMappedOffsetInletFvPatchScalarField_H
#define monitorMappedOffsetInletFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class monitorMappedOffsetInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class monitorMappedOffsetInletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Name of the patch to be monitored
        word monitorPatchName_;

        //- Name of flux field
        word phiName_;

        //- Scalar value to add or subtract
        scalar offset_;

        //- Upper-limit scalar
        scalar max_;

        //- Lower-limit scalar
        scalar min_;


public:

    //- Runtime type information
    TypeName("monitorMappedOffsetInlet");


    // Constructors

        //- Construct from patch and internal field
        monitorMappedOffsetInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        monitorMappedOffsetInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- monitorMappedOffsetInletFvPatchScalarField
        //- onto a new patch
        monitorMappedOffsetInletFvPatchScalarField
        (
            const monitorMappedOffsetInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        monitorMappedOffsetInletFvPatchScalarField
        (
            const monitorMappedOffsetInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new monitorMappedOffsetInletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        monitorMappedOffsetInletFvPatchScalarField
        (
            const monitorMappedOffsetInletFvPatchScalarField&,
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
                new monitorMappedOffsetInletFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

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
