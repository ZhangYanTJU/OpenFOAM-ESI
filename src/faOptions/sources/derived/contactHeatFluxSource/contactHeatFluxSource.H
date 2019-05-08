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
    Foam::fa::contactHeatFluxSource

Group
    grpFvOptionsSources

Description
    Contact heat flux between faMesh and fvMesh


Usage


SourceFiles
    contactHeatFluxSource.C


\*---------------------------------------------------------------------------*/

#ifndef fa_contactHeatFluxSource_H
#define fa_contactHeatFluxSource_H

#include "faOption.H"
#include "Function1.H"
#include "areaFields.H"
#include "faceSetOption.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                        Class contactHeatFluxSource Declaration
\*---------------------------------------------------------------------------*/

class contactHeatFluxSource
:
    public faceSetOption,
    public temperatureCoupledBase
{
    // Private data

        //- Name of temperature field - default = "T" (optional)
        word TName_;

        //- Name of primary temperature field
        word TprimaryName_;

        //- Primary region temperature
        const volScalarField& Tp_;

        //- Temperature - wall [K]
        areaScalarField Tw1_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;

        //- Total contact resistance
        scalar contactRes_;

        //- Current time index (used for updating)
        label curTimeIndex_;


    // Private Member Functions

        //- No copy construct
        contactHeatFluxSource(const contactHeatFluxSource&) = delete;

        //- No copy assignment
        void operator=(const contactHeatFluxSource&) = delete;

        //- Return htc from the primary region
        tmp<areaScalarField> htc() const;


public:

    //- Runtime type information
    TypeName("contactHeatFluxSource");


    // Constructors

        //- Construct from explicit source name and mesh
        contactHeatFluxSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const fvPatch& patch
        );


    //- Destructor
    virtual ~contactHeatFluxSource();


    // Member Functions

        // Evaluate

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const areaScalarField& h,
                const areaScalarField& rho,
                faMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
