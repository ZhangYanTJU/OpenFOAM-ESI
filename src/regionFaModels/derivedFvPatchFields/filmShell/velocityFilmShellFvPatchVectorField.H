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
    Foam::compressible::velocityFilmShellFvPatchVectorField

Group
    grpLiquidFilmBoundaryConditions

Description


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
        type                velocityFilmShell;

        active              true;

        U                   U;

        thermo
        {

        }

        region              ceilingShell;
        liquidFilm          kinematicThinFilm;

        value               uniform (0 0 0);

    }
    \endverbatim


SourceFiles
    velocityFilmShellFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef velocityFilmShellFvPatchVectorField_H
#define velocityFilmShellFvPatchVectorField_H


#include "autoPtr.H"
#include "liquidFilmBase.H"
#include "mixedFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class velocityFilmShellFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class velocityFilmShellFvPatchVectorField
:
    public mixedFvPatchField<vector>
{

    // Private data

        //- Thermal baffle
        autoPtr<regionModels::areaSurfaceFilmModels::liquidFilmBase> baffle_;

        //- Dictionary
        mutable dictionary dict_;

        //- Time index to evolve the film
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("velocityFilmShell");


    // Constructors

        //- Construct from patch and internal field
        velocityFilmShellFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        velocityFilmShellFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  velocityFilmShellFvPatchVectorField onto a new patch
        velocityFilmShellFvPatchVectorField
        (
            const velocityFilmShellFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new velocityFilmShellFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        velocityFilmShellFvPatchVectorField
        (
            const velocityFilmShellFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new velocityFilmShellFvPatchVectorField(*this, iF)
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
