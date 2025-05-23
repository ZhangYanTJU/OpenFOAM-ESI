/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::externalWallHeatFluxTemperatureFvPatchScalarField

Group
    grpThermoBoundaryConditions grpWallBoundaryConditions

Description
    This boundary condition applies a heat flux condition to temperature
    on an external wall in one of three modes:

      - fixed power: supply Q
      - fixed heat flux: supply q
      - fixed heat transfer coefficient: supply h and Ta

    where:
    \vartable
        Q  | Power [W]
        q  | Heat flux [W/m^2]
        h  | Heat transfer coefficient [W/m^2/K]
        Ta | Ambient temperature [K]
    \endvartable

    For heat transfer coefficient mode optional thin thermal layer resistances
    can be specified through thicknessLayers and kappaLayers entries.

    The thermal conductivity \c kappa can either be retrieved from various
    possible sources, as detailed in the class temperatureCoupledBase.

    The ambient temperature Ta is specified as a Foam::Function1 of time but
    uniform in space.

Usage
    \table
    Property    | Description                              | Required | Default
    mode        | 'power', 'flux' or 'coefficient'         | yes |
    Q           | Power [W]                                | for mode 'power' |
    q           | Heat flux [W/m^2]                        | for mode 'flux' |
    h           | Heat transfer coefficient [W/m^2/K] | for mode 'coefficient' |
    Ta          | Ambient temperature [K]             | for mode 'coefficient' |
    thicknessLayers | Layer thicknesses [m]                | no  |
    kappaLayers | Layer thermal conductivities [W/m/K]     | no  |
    relaxation  | Relaxation for the wall temperature      | no  | 1
    emissivity  | Surface emissivity for radiative flux to ambient | no | 0
    qr          | Name of the radiative field              | no  | none
    qrRelaxation | Relaxation factor for radiative field   | no  | 1
    kappaMethod | Inherited from temperatureCoupledBase    | inherited |
    kappa       | Inherited from temperatureCoupledBase    | inherited |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            externalWallHeatFluxTemperature;

        mode            coefficient;

        Ta              constant 300.0;
        h               constant 10.0;
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);

        kappaMethod     fluidThermo;

        value           $internalField;
    }
    \endverbatim

Note
    Quantities that are considered "global" (eg, power, ambient temperature)
    can be specified as Function1 types.
    Quantities that may have local variations (eg, htc, heat-flux)
    can be specified as PatchFunction1 types.

See also
    Foam::temperatureCoupledBase
    Foam::mixedFvPatchScalarField

SourceFiles
    externalWallHeatFluxTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef externalWallHeatFluxTemperatureFvPatchScalarField_H
#define externalWallHeatFluxTemperatureFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class externalWallHeatFluxTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class externalWallHeatFluxTemperatureFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    // Public Data

        //- Operation mode enumeration
        enum operationMode
        {
            fixedPower,             //!< Heat power [W]
            fixedHeatFlux,          //!< Heat flux [W/m2]
            fixedHeatTransferCoeff, //!< Heat transfer coefficient [W/m^2/K]
        };

        static const Enum<operationMode> operationModeNames;


private:

    // Private Data

        //- Operation mode
        operationMode mode_;

        //- Heat power [W]
        autoPtr<Function1<scalar>> Q_;

        //- Heat flux [W/m2]
        autoPtr<PatchFunction1<scalar>> q_;

        //- Heat transfer coefficient [W/m2K]
        autoPtr<PatchFunction1<scalar>> h_;

        //- Ambient temperature [K]
        autoPtr<Function1<scalar>> Ta_;

        //- Relaxation for the wall temperature (thermal inertia)
        scalar relaxation_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar emissivity_;

        //- Cache qr for relaxation
        scalarField qrPrevious_;

        //- Relaxation for qr
        scalar qrRelaxation_;

        //- Name of the radiative heat flux
        const word qrName_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;


public:

    //- Runtime type information
    TypeName("externalWallHeatFluxTemperature");


    // Constructors

        //- Construct from patch and internal field
        externalWallHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        externalWallHeatFluxTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  externalWallHeatFluxTemperatureFvPatchScalarField
        //  onto a new patch
        externalWallHeatFluxTemperatureFvPatchScalarField
        (
            const externalWallHeatFluxTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        externalWallHeatFluxTemperatureFvPatchScalarField
        (
            const externalWallHeatFluxTemperatureFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        externalWallHeatFluxTemperatureFvPatchScalarField
        (
            const externalWallHeatFluxTemperatureFvPatchScalarField&,
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

        // Access

            //- Allow manipulation of the boundary values
            virtual bool fixesValue() const
            {
                return false;
            }


        // Mapping functions

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
