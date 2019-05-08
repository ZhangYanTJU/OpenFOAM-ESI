/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Foam::compressible::thermalShellFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a coupled temperature condition between
    a primary region (3D mesh) and a thermal shell model (2D mesh).


    The primary region boundary creates the finite are region and evolves
    its energy equation.

    NOTE: The 2D area mesh needs to be generate in the pre-processing steps


Usage
    \table
    Property    | Description                              | Required | Default
    kappaMethod | 'fluidThermo', 'solidThermo'             | yes
    thermo      | Solid thermo                             | yes
    thermalShellModel  | thermal shell modes               | yes
    region      | 2D region to be applied                  | yes
    T           | Temperature field name                   | yes
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <masterPatchName>
    {
        type                compressible::thermalShell;
        kappaMethod         fluidThermo;
        kappa               none;

        active              true;

        T                   T;

        thermo
        {
            rho                 1000;
            kappa               200;
            Cp                  600;
            Hf                  0;
            emissivity          0;
        }

        region              ceilingShell;
        thermalShellModel   thermalShell;

        value           uniform 300;

    }
    \endverbatim


SourceFiles
    thermalShellFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermalShellFvPatchScalarField_H
#define thermalShellFvPatchScalarField_H


#include "autoPtr.H"
#include "thermalShellModel.H"
#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
             Class thermalShellFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class thermalShellFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
    // Private data

        //- Thermal baffle
        autoPtr<regionModels::thermalShellModel> baffle_;

        //- Dictionary
        mutable dictionary dict_;


public:

    //- Runtime type information
    TypeName("compressible::thermalShell");


    // Constructors

        //- Construct from patch and internal field
        thermalShellFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalShellFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  thermalShellFvPatchScalarField onto a new patch
        thermalShellFvPatchScalarField
        (
            const thermalShellFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalShellFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thermalShellFvPatchScalarField
        (
            const thermalShellFvPatchScalarField&,
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
                new thermalShellFvPatchScalarField(*this, iF)
            );
        }


    // Member functions


        // Mapping functions
/*
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

*/
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
