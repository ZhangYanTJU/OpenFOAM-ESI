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
    Foam::slowSorptionFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a mixed-value condition for
    a given scalar field to model kinetically-controlled slow sorption
    processes which PENDING (P:p. 18).
    Available slow sorption models are elaborated below.

    Boundary condition expression (C:Eq. 3.21):
    \f[
        D \vec{n} \cdot (\grad c)_p = -S_p
    \f]

    where
    \vartable
      \D   | Diffusion constant
      c_p  | Adsorbed substance profile on patch
      c_o  | Substance profile in bulk flow
      S_p  | Net rate of reversible adsorption on patch
    \endvartable

    Henry slow sorption model (P:Eq. 2.55):

    \f[
        S_p = \kappa^{ads} c_{sub} - \kappa^{des} \frac{c_p}{c_{p_{max}}}
    \f]

    where
    \vartable
      \kappa^{ads} | Adsorption rate                      [m/s]
      \kappa^{des} | Desorption rate                      [x/(m^2 s)]
      c_p          | Adsorbed substance profile on patch
      c_{p_{max}}  | Maximum substance that can be adsorbed on patch
      c_{sub}      | Substance profile on subsurface between bulk flow and patch
    \endvartable


    Langmuir slow sorption model (P:Eq. 2.56):

    \f[
        S_p =
            \kappa^{ads} c_{sub} \left(1 - \frac{c_p}{c_{p_{max}}} \right)
          - \kappa^{des} \frac{c_p}{c_{p_{max}}}
    \f]

    References:
    \verbatim
        Sorption boundary condition and slow-sorption models (tag:P):
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
        type        slowSorption;
        model       <slowSorptionModel>;

        D           <scalar>;
        kAds        <scalar>;
        max         <scalar>;

        // Optional entries
        kDes        <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                      | Type   | Reqd   | Deflt
      type     | Type name: slowSorption          | word   | yes    | -
      model    | Slow sorption model              | word   | yes    | -
      D        | Diffusion coefficient            | scalar | yes    | -
      kAds     | Adsorption rate                  | scalar | yes    | -
      max      | Maximum level of adsorption of a given substance on patch <!--
               -->                                | scalar | cndtnl | -
      kDes     | Desorption rate                  | scalar | no     | 0
    \endtable

    The inherited entries are elaborated in:
      - \link mixedFvPatchFields.H \endlink

    Options for the \c model entry:
    \verbatim
      Henry            | Henry slow sorption model
      Langmuir         | Langmuir slow sorption model
    \endverbatim

Note
  - The subsurface concentration is presumed to be equal to the bulk
    concentration by assuming the adsorption rate (\c kAds) is from
    bulk flow to patch rather than the bulk flow to subsurface (P:Fig. 2.2).

SourceFiles
    slowSorptionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef slowSorptionFvPatchScalarField_H
#define slowSorptionFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class slowSorptionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class slowSorptionFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Enumerations

        //- Options for the slow sorption model
        enum slowSorptionModelType : char
        {
            HENRY = 0,          //!< "Henry fast-sorption model"
            LANGMUIR            //!< "Langmuir fast-sorption model"
        };

        //- Names for slowSorptionModelType
        static const Enum<slowSorptionModelType> slowSorptionModelTypeNames;


    // Private Data

        //- Sorption model
        enum slowSorptionModelType slowSorptionModel_;

        //- Diffusion coefficient
        scalar D_;

        //- Adsorption rate
        scalar kAds_;

        //- Maximum level of adsorption of a given substance on patch
        scalar max_;

        //- Desorption rate
        scalar kDes_;


public:

    //- Runtime type information
    TypeName("slowSorption");


    // Constructors

        //- Construct from patch and internal field
        slowSorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        slowSorptionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- slowSorptionFvPatchScalarField onto a new patch
        slowSorptionFvPatchScalarField
        (
            const slowSorptionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        slowSorptionFvPatchScalarField
        (
            const slowSorptionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new slowSorptionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        slowSorptionFvPatchScalarField
        (
            const slowSorptionFvPatchScalarField&,
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
                new slowSorptionFvPatchScalarField(*this, iF)
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
