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
    Foam::meltingEvaporationModels::interfaceOxide

Description



Usage

    Example usage:
    \verbatim
    massTransferModel
    (
        (liquid to gas)
        {
            type                interfaceOxide;

            C                   0.1;

            Tactivate           373;
        }
    );
    \endverbatim

    where:
    \table
        Property    | Description             | Required    | Default value
        C           | Coefficient (C > 0 for evaporation, C < 0 for
            condensation) | yes
        includeVolChange    | Volumen change  | no          | yes
        isoAlpha    | iso-alpha for interface | no          | 0.5
        Tactivate   | Saturation temperature  | yes
        species     | Specie name on the other phase | no   | none
    \endtable

SourceFiles
    interfaceOxide.C

\*---------------------------------------------------------------------------*/

#ifndef meltingEvaporationModels_interfaceOxide_H
#define meltingEvaporationModels_interfaceOxide_H

#include "InterfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{

// Forward Declarations
class phasePair;

namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                    Class interfaceOxide Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class interfaceOxide
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private Data

        //- Model constant
        dimensionedScalar C_;

        //- Activation temperature
        const dimensionedScalar Tactivate_;

        //- Critical alpha
        dimensionedScalar alphaCrit_;

        //- Interface area
        //volScalarField interfaceArea_;

        //- Heat transfer coefficient
        //volScalarField htc_;

        //- Mass source
        volScalarField mDotc_;

        //- Interface Iso-value
        //scalar isoAlpha_;


    // Private Member Functions

        //- Update interface
        //void updateInterface(const volScalarField& T);


public:

    //- Runtime type information
    TypeName("interfaceOxideRate");


    // Constructors

        //- Construct from components
        interfaceOxide
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interfaceOxide() = default;


    // Member Functions

        //- Explicit total mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            const volScalarField& field
        );

        //- Implicit mass transfer coefficient
        virtual tmp<volScalarField> KSp
        (
            label modelVariable,
            const volScalarField& field
        );

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> KSu
        (
            label modelVariable,
            const volScalarField& field
        );

        //- Return Tactivate
        virtual const dimensionedScalar& Tactivate() const;

        //- Add/subtract alpha*div(U) as a source term
        //- for alpha, substituting div(U) = mDot(1/rho1 - 1/rho2)
        virtual bool includeDivU();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interfaceOxide.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
