/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::enthalpySorptionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This is a temperature BC which works in conjunction with speciesSorption
    for a species.

    This boundary condition substract or adds enthalpy associated with the
    adsorption provided by the BC speciesSorption.

    It can handle two enthalpy models:

    1) Estimate
    2) Calculated

    On top of this the enthalphy associated with the sensible enthalphy
    corresponding with the species transport can be added using includeHs


Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type                enthalpySorption;
        enthalpyModel       calculated; //estimated;
        enthalpyTable
        {
            type            table;
            values          ((0 0)(1 50));
        }
        C                   2;
        species             O2;
        includeHs           true;
        value               $internalField;
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                        | Type   | Reqd   | Deflt
      type     | Type name: enthalpySorption        | word   | yes    | -
      C        | Estimate enthalphy model  constant | scalar | yes    | -
      species  | Species associated with this BC    | word   | yes    | -
      includeHs| Include sensible enthalpy          | bool   | no     | true
      enthalpyModel| Adsorption enthalphy model     | word   | yes    | -
      enthalpyTable| Calculated enthalphy model table | Function1 | no
      Hvap         | Evaporation heat for species   | scalar | no

    \endtable


SourceFiles
    enthalpySorptionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef enthalpySorptionFvPatchScalarField_H
#define enthalpySorptionFvPatchScalarField_H

#include "boundarySourcePatch.H"
#include "zeroGradientFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class enthalpySorptionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class enthalpySorptionFvPatchScalarField
:
    public zeroGradientFvPatchField<scalar>,
    public boundarySourcePatch
{
private:

    // Private Enumeration


        //- Options for the enthalpy model
        enum enthalpyModelType
        {
            estimated,
            calculated
        };


        //- Names for kineticModelType
        static const Enum<enthalpyModelType> enthalpyModelTypeNames;


private:

    // Private Data


        //- Enthalpy model
        enum enthalpyModelType enthalpyModel_;

        //- Estimated enthalpy model constant
        scalar C_;

        //- Load enthalpy Table for calculated model
        autoPtr<Function1<scalar>> enthalpyMassLoad_;

        //- Species field name
        word speciesName_;

        //- Enthalpy change on cells next to patch [J/Kg]
        scalarField dhdt_;

        //- Include sensible enthalpy of the species
        bool includeHs_;

        //- Pressure field name
        word pName_;

        //- Temperature field name
        word TName_;

        //- Heat of evaporation of species
        scalar Hvap_;


public:

    //- Runtime type information
    TypeName("enthalpySorption");


    // Constructors

        //- Construct from patch and internal field
        enthalpySorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        enthalpySorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- enthalpySorptionFvPatchScalarField onto a new patch
        enthalpySorptionFvPatchScalarField
        (
            const enthalpySorptionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        enthalpySorptionFvPatchScalarField
        (
            const enthalpySorptionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new enthalpySorptionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        enthalpySorptionFvPatchScalarField
        (
            const enthalpySorptionFvPatchScalarField&,
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
                new enthalpySorptionFvPatchScalarField(*this, iF)
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
