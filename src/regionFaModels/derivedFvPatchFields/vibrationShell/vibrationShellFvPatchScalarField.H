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
    Foam::compressible::vibrationShellFvPatchScalarField

Group
    grpVibrationBoundaryConditions

Description
   


Usage
    \table
    Property    | Description                              | Required | Default
    thermo      | Solid thermo                             | yes
    vibrationShellModel  | vibration shell model           | yes
    region      | 2D region to be applied                  | yes
    p           | Pressure field name                      | yes
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <masterPatchName>
    {
        type                vibrationShell;
        active              true;
        p                   pa;

        thermo
        {
            rho                 1000;
            kappa               200;
            Cp                  600;
            Hf                  0;
            emissivity          0;
            
            E                   0;
            nu                  1;
        }

        region                ceilingShell;
        vibrationShellModel   KirchhoffShell;

        value           uniform 0;

    }
    \endverbatim


SourceFiles
    vibrationShellFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef vibrationShellFvPatchScalarField_H
#define vibrationShellFvPatchScalarField_H


#include "autoPtr.H"
#include "vibrationShellModel.H"
#include "mixedFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class vibrationShellFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class vibrationShellFvPatchScalarField
:
    public mixedFvPatchField<scalar>
{
    // Private data

        //- Thermal baffle
        autoPtr<regionModels::vibrationShellModel> baffle_;

        //- Dictionary
        mutable dictionary dict_;


public:

    //- Runtime type information
    TypeName("vibrationShell");


    // Constructors

        //- Construct from patch and internal field
        vibrationShellFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        vibrationShellFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  vibrationShellFvPatchScalarField onto a new patch
        vibrationShellFvPatchScalarField
        (
            const vibrationShellFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new vibrationShellFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        vibrationShellFvPatchScalarField
        (
            const vibrationShellFvPatchScalarField&,
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
                new vibrationShellFvPatchScalarField(*this, iF)
            );
        }


    // Member functions


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
