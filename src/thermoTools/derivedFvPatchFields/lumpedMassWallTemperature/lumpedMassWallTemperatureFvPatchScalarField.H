/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::lumpedMassWallTemperatureFvPatchScalarField

Group
    grpThermoBoundaryConditions grpWallBoundaryConditions

Description
    Employs a lumped mass model for temperature.

    It considers a single temperature value for the whole patch and evaluates
    the temperature evolution using the net heat flux into the patch.

Usage
    \table
        Property     | Description                 | Required | Default value
        kappaMethod  | thermal conductivity option | yes      |
        Cp           | Heat capacity [J/kg.K]      | yes      |
        mass         | Total mass [kg]             | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            lumpedMassWallTemperature;
        kappaMethod     fluidThermo;
        kappa           none;
        mass            1000;
        Cp              4100;
        value           uniform 300.0;
    }
    \endverbatim

SourceFiles
    lumpedMassWallTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedMassWallTemperatureFvPatchScalarField_H
#define lumpedMassWallTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class lumpedMassWallTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class lumpedMassWallTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
private:

    // Private data

        //- Specific heat capacity [J/kg/K]
        scalar Cp_;

        //- Mass [kg]
        scalar mass_;

        //- Current time index (used for updating)
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("lumpedMassWallTemperature");


    // Constructors

        //- Construct from patch and internal field
        lumpedMassWallTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        lumpedMassWallTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // lumpedMassWallTemperatureFvPatchScalarField
        // onto a new patch
        lumpedMassWallTemperatureFvPatchScalarField
        (
            const lumpedMassWallTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        lumpedMassWallTemperatureFvPatchScalarField
        (
            const lumpedMassWallTemperatureFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        lumpedMassWallTemperatureFvPatchScalarField
        (
            const lumpedMassWallTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
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
                const fvPatchField<scalar>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
