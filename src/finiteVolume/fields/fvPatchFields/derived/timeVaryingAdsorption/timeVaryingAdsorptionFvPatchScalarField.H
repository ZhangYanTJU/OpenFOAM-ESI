/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::timeVaryingAdsorptionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a time-varying fixed-value condition
    for a given scalar field to model time-dependent adsorption processes:

    \f[
        \frac{d c}{d t} = k (c_{p_{max}} - c)
    \f]

    where
    \vartable
      c            | Adsorbed substance on patch at time t
      k            | Adsorption rate constant
      c_{p_{max}}  | Maximum substance that can be adsorbed on patch
    \endvartable

    References:
    \verbatim
        Pseudo-first order adsorption model (tag:ANK):
            Adelodun, A. A., Ngila, J. C., Kim, D. G., & Jo, Y. M. (2016).
            Isotherm, thermodynamic and kinetic studies of selective
            CO2 adsorption on chemically modified carbon surfaces.
            Aerosol and Air Quality Research, 16(12), 3312-3329.
            DOI:10.4209/aaqr.2016.01.0014
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type        timeVaryingAdsorption;
        k           <scalar>;
        max         <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                      | Type   | Reqd   | Deflt
      type     | Type name: timeVaryingAdsorption | word   | yes    | -
      k        | Adsorption rate constant         | scalar | yes    | -
      max      | Maximum level of adsorption of a given substance on patch <!--
               -->                                | scalar | yes    | -
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink

SourceFiles
    timeVaryingAdsorptionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingAdsorptionFvPatchScalarField_H
#define timeVaryingAdsorptionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class timeVaryingAdsorptionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingAdsorptionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
public:

    // Public Enumeration

        //- Enumeration defining the available ddt schemes
        enum ddtSchemeType
        {
            tsEuler,
            tsCrankNicolson,
            tsBackward
        };


private:

    // Private Data

        //- Time scheme type names
        static const Enum<ddtSchemeType> ddtSchemeTypeNames_;

        //- Adsorption rate constant
        scalar k_;

        //- Maximum level of adsorption of a given substance on patch
        scalar max_;


public:

    //- Runtime type information
    TypeName("timeVaryingAdsorption");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingAdsorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingAdsorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- timeVaryingAdsorptionFvPatchScalarField onto a new patch
        timeVaryingAdsorptionFvPatchScalarField
        (
            const timeVaryingAdsorptionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingAdsorptionFvPatchScalarField
        (
            const timeVaryingAdsorptionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new timeVaryingAdsorptionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingAdsorptionFvPatchScalarField
        (
            const timeVaryingAdsorptionFvPatchScalarField&,
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
                new timeVaryingAdsorptionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation

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
