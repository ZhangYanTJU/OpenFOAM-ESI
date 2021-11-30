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
    Foam::speciesSorptionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a first order zero gradient condition
    for a given scalar field to model time-dependent adsorption-desoprtion

    \f[
        \frac{d c}{d t} = kasb (c_{eq} - c_{abs})
    \f]


    where
    \vartable
      c_{eq}       | Equilibrium concentration
      c_{abs}      | Absorbed at wall
      kasb         | Adsorption rate constant [1/s]
    \endvartable

    \f[
        c_{eq} = max* (kl_*co_{int}/(1 + kl_*co_{int}));
    \f]

    where
    \vartable
      max          | MAximum concentration
      kl           | Langmuir constant
      co_{int}     | Local cell value concentration
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type        speciesSorption;
        kbas        <scalar>;
        max         <scalar>;
        model       Langmuir;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                        | Type   | Reqd   | Deflt
      type     | Type name: speciesSorption | word  | yes    | -
      kbas     | Adsorption rate constant           | scalar | yes    | -
      kl       | Langmuir constant                  | scalar | yes    | -
      max      | Maximum concentation at wall       | scalar | yes    | -
      model    | Equilibrium model                  | word   | yes
    \endtable


SourceFiles
    speciesSorptionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef speciesSorptionFvPatchScalarField_H
#define speciesSorptionFvPatchScalarField_H

#include "boundarySourcePatch.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class speciesSorptionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class speciesSorptionFvPatchScalarField
:
    public zeroGradientFvPatchField<scalar>,
    public boundarySourcePatch
{
public:

    // Public Enumeration

        //- Options for the slow sorption model
        enum equilibriumModelType : char
        {
            LANGMUIR = 0
        };


        //- Names for slowSorptionModelType
        static const Enum<equilibriumModelType> equilibriumModelTypeNames;


private:

    // Private Data

         //- Equilibrium model
        enum equilibriumModelType equilibriumModel_;

        //- Adsorption rate constant  [1/sec]
        scalar kabs_;

        //- Langmuir adsorption constant
        scalar kl_;

        //- Maximum level of adsorption of a given substance on patch [mol/Kg]
        scalar max_;

        //- Source on cells next to patch [mol/Kg/sec]
        scalarField dfldp_;

        //- Absorbed mass
        scalarField mass_;

    // Private Member Functions

        //- Calculate the mole fraction fields
        tmp<scalarField> calcMoleFractions();


public:

    //- Runtime type information
    TypeName("speciesSorption");


    // Constructors

        //- Construct from patch and internal field
        speciesSorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        speciesSorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- speciesSorptionFvPatchScalarField onto a new patch
        speciesSorptionFvPatchScalarField
        (
            const speciesSorptionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        speciesSorptionFvPatchScalarField
        (
            const speciesSorptionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new speciesSorptionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        speciesSorptionFvPatchScalarField
        (
            const speciesSorptionFvPatchScalarField&,
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
                new speciesSorptionFvPatchScalarField(*this, iF)
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

            //- Source of cells next to the patch
            virtual tmp<scalarField> patchSource() const;


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
