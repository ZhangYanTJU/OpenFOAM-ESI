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
    Foam::fastSorptionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a fixed-value condition for
    a given scalar field to model diffusion-controlled fast sorption
    processes which are characterised by quasi-instantaneous sorption
    (P:p. 18). Available fast sorption models are elaborated below.

    Henry fast sorption model (P:Eq. 2.39, p. 18), with which a linear
    relationship between bulk concentration \f$ c_o \f$ and concentration
    of substance adsorbed on patch \f$ c_p \f$ is assumed:

    \f[
        c_p = k c_o
    \f]

    where
    \vartable
      c_p      | Adsorbed substance profile on patch
      k        | Henry constant
      c_o      | Substance profile in bulk flow
    \endvartable


    Langmuir fast sorption model (P:Eq. 2.41, p. 18), with which
    interactions on patch are partially considered:

    \f[
        c_p = c_{p_{max}} \frac{c_o}{a_L + c_o}
    \f]

    where
    \vartable
      c_{p_{max}}    | Maximum substance that can be adsorbed on patch
      a_L            | Langmuir equilibrium constant
    \endvartable

    References:
    \verbatim
        Fast sorption models (tag:P):
            Pesci, C. (2019).
            Computational analysis of fluid interfaces
            influenced by soluble surfactant.
            PhD thesis. Darmstadt: Technische Universitat Darmstadt.
            URI:https://tuprints.ulb.tu-darmstadt.de/id/eprint/9303
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type        fastSorption;
        model       <fastSorptionModel>;

        // Conditional mandatory entries

            // Option-1
            HenryCoeffs
            {
                k           <scalar>;
            }

            // Option-2
            LangmuirCoeffs
            {
                max         <scalar>;
                aL          <scalar>;
            }

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                      | Type   | Reqd   | Deflt
      type     | Type name: fastSorption          | word   | yes    | -
      model    | Fast sorption model              | word   | yes    | -
      k        | Henry constant                   | scalar | cndtnl | -
      max      | Maximum level of adsorption of a given substance on patch <!--
               -->                                | scalar | cndtnl | -
      aL       | Langmuir equilibrium constant    | scalar | cndtnl | -
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink

    Options for the \c model entry:
    \verbatim
      Henry            | Henry fast sorption model
      Langmuir         | Langmuir fast sorption model
    \endverbatim

SourceFiles
    fastSorptionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fastSorptionFvPatchScalarField_H
#define fastSorptionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class fastSorptionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fastSorptionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Enumerations

        //- Options for the fast sorption model
        enum fastSorptionModelType : char
        {
            HENRY = 0,          //!< "Henry fast-sorption model"
            LANGMUIR            //!< "Langmuir fast-sorption model"
        };

        //- Names for fastSorptionModelType
        static const Enum<fastSorptionModelType> fastSorptionModelTypeNames;


    // Private Data

        //- Sorption model
        enum fastSorptionModelType fastSorptionModel_;


        // Henry fast-sorption model coefficients

            //- Henry constant
            scalar k_;


        // Langmuir fast-sorption model coefficients

            //- Maximum level of adsorption of a given substance on patch
            scalar max_;

            //- Langmuir equilibrium constant
            scalar aL_;


public:

    //- Runtime type information
    TypeName("fastSorption");


    // Constructors

        //- Construct from patch and internal field
        fastSorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fastSorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- fastSorptionFvPatchScalarField onto a new patch
        fastSorptionFvPatchScalarField
        (
            const fastSorptionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fastSorptionFvPatchScalarField
        (
            const fastSorptionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fastSorptionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fastSorptionFvPatchScalarField
        (
            const fastSorptionFvPatchScalarField&,
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
                new fastSorptionFvPatchScalarField(*this, iF)
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
