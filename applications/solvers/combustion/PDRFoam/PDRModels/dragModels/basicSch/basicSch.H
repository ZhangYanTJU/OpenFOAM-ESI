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
    Foam::PDRDragModels::basicSch

Description


SourceFiles
    basicSch.C

\*---------------------------------------------------------------------------*/

#ifndef basicSch_H
#define basicSch_H

#include "PDRDragModel.H"
#include "XiEqModel.H"


namespace Foam
{
namespace PDRDragModels
{

/*---------------------------------------------------------------------------*\
                          Class basicSch Declaration
\*---------------------------------------------------------------------------*/

class basicSch
:
    public PDRDragModel
{
    // Private data

        dimensionedScalar Csu;
        dimensionedScalar Csk;
        volScalarField Aw_;
        volSymmTensorField CR_;

        //- Schelkin effect Model constants
        const scalar nrCoef_;
        const scalar nrExp2_;
        const scalar lCoef_;
        const scalar maxSchFac_;

        //- Use sub-grid Schelkin effect
        bool subGridSchelkin_;


    // Private Member Functions

        //- Disallow copy construct
        basicSch(const basicSch&);

        //- Disallow default bitwise assignment
        void operator=(const basicSch&);


public:

    //- Runtime type information
    TypeName("basicSch");


    // Constructors

        //- Construct from components
        basicSch
        (
            const dictionary& PDRProperties,
            const compressible::RASModel& turbulence,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~basicSch();


    // Member Functions

        //- Return the momentum drag coefficient
        virtual tmp<volSymmTensorField> Dcu() const;

        //- Return the momentum drag turbulence generation rate
        virtual tmp<volScalarField> Gk() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& PDRProperties);

        //- Write fields
        void writeFields() const;

        //- Return the schelkin factor for drag turbulence generation rate
        tmp<volScalarField> schFac() const;

        //- Return the sub-grid Schelkin effect exponent
        tmp<volScalarField> SchelkinExponent
        (
            const scalar,
            const scalar,
            const volScalarField&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PDRDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
