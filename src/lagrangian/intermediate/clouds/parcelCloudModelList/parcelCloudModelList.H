/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#ifndef Foam_parcelCloudModelList_H
#define Foam_parcelCloudModelList_H

#include "parcelCloudModel.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class parcelCloudModelList Declaration
\*---------------------------------------------------------------------------*/

class parcelCloudModelList
:
    public PtrList<parcelCloudModel>
{
private:

    //- Reference to the mesh
    const fvMesh& mesh_;


public:

    // Constructors

        parcelCloudModelList
        (
            const dimensionedVector& g,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu
        );


        parcelCloudModelList
        (
            const dimensionedVector& g,
            const volScalarField& rho,
            const volVectorField& U,
            const SLGThermo& slgThermo
        );


    // Member Functions

        // Evolution

            //- Evolve the cloud collection
            void evolve();


        // Source terms

            //- Return momentum transfer [kg m/s]
            inline tmp<volVectorField::Internal> UTrans() const;

            //- Return momentum source term [kg m/s^2]
            inline tmp<fvVectorMatrix> SU(volVectorField& U) const;


            //- Return sensible enthalpy transfer [J]
            inline tmp<volScalarField::Internal> hsTrans() const;

            //- Return sensible enthalpy source term [J/s]
            inline tmp<fvScalarMatrix> Sh(volScalarField& hs) const;


            //- Return mass source term for specie [kg/s]
            inline tmp<fvScalarMatrix> SYi
            (
                const label speciei,
                volScalarField& Yi
            ) const;

            //- Return total mass source term [kg/m^3/s]
            inline tmp<volScalarField::Internal> Srho() const;

            //- Return total mass source term [kg/m3/s]
            inline tmp<fvScalarMatrix> Srho(volScalarField& rho) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "parcelCloudModelListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
