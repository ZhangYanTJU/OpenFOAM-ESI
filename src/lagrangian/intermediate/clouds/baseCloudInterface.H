/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::baseCloudInterface

Description
    Virtual abstract base class for parcel clouds. Inserted by ParcelCloudModel
    into the base of the cloud template hierarchy and adds virtualisation of
    most methods defined by the clouds.

    Note that the "evolve" method is not virtualised here. Due to the way in
    which TrackCloudType and trackingData templating work, it is not possible
    to virtualise this method directly. Instead it has to be wrapped. That is
    achieved by the parcelCloudModel and ParcelCloudModel classes.

SourceFiles
    baseCloudInterface.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_baseCloudInterface_H
#define Foam_baseCloudInterface_H

#include "volFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class baseCloudInterface
{
public:

    //- Runtime type information
    TypeName("baseCloudInterface");

    // Null constructor
    baseCloudInterface()
    {}

    //- Destructor
    virtual ~baseCloudInterface() = default;


    // Member Functions

        //- Number of parcels
        virtual label nParcels() const = 0;

        //- Total mass in system
        virtual scalar massInSystem() const = 0;

        //- Total linear momentum of the system
        virtual vector linearMomentumOfSystem() const = 0;

        //- Total linear kinetic energy in the system
        virtual scalar linearKineticEnergyOfSystem() const = 0;

        //- Mean diameter Dij
        virtual scalar Dij(const label i, const label j) const = 0;

        //- Max diameter
        virtual scalar Dmax() const = 0;


        // Fields

            //- Volume swept rate of parcels per cell
            virtual const tmp<volScalarField> vDotSweep() const = 0;

            //- Return particle volume fraction field
            //  Note: for particles belonging to this cloud only
            virtual const tmp<volScalarField> theta() const = 0;

            //- Return particle mass fraction field
            //  Note: for particles belonging to this cloud only
            virtual const tmp<volScalarField> alpha() const = 0;

            //- Return particle effective density field
            //  Note: for particles belonging to this cloud only
            virtual const tmp<volScalarField> rhoEff() const = 0;


        // Sources

            // Momentum

                // //- Return reference to momentum source
                // virtual volVectorField::Internal& UTrans() = 0;

                // //- Return const reference to momentum source
                // virtual const volVectorField::Internal& UTrans() const = 0;

                // //- Return coefficient for carrier phase U equation
                // virtual volScalarField::Internal& UCoeff() = 0;

                // //- Return const coefficient for carrier phase U equation
                // virtual const volScalarField::Internal& UCoeff() const = 0;

                //- Return momentum source term (compressible)
                virtual tmp<fvVectorMatrix> SU(volVectorField& U) const = 0;


            // Energy

                // //- Sensible enthalpy transfer [J/kg]
                // virtual volScalarField::Internal& hsTrans() = 0;

                // //- Sensible enthalpy transfer [J/kg]
                // virtual const volScalarField::Internal&
                //     hsTrans() const = 0;

                // //- Return coefficient for carrier phase hs equation
                // virtual volScalarField::Internal& hsCoeff() = 0;

                // //- Return const coefficient for carrier phase hs equation
                // virtual const volScalarField::Internal&
                //     hsCoeff() const = 0;

                //- Return sensible enthalpy source term [J/s]
                virtual tmp<fvScalarMatrix> Sh(volScalarField& hs) const = 0;


                // Radiation - implemented by ThermoCloud level

                    //- Return equivalent particulate emission
                    virtual tmp<volScalarField> Ep() const = 0;

                    //- Return equivalent particulate absorption
                    virtual tmp<volScalarField> ap() const = 0;

                    //- Return equivalent particulate scattering factor
                    virtual tmp<volScalarField> sigmap() const = 0;


            // Mass

                //- Return mass source term for specie i - specie eqn
                virtual tmp<fvScalarMatrix> SYi
                (
                    const label i,
                    volScalarField& Yi
                ) const = 0;

                //- Return mass source for field i - fully explicit
                virtual tmp<volScalarField::Internal> Srho
                (
                    const label i
                ) const = 0;

                //- Return total mass source for carrier phase
                virtual tmp<volScalarField::Internal> Srho() const = 0;

                //- Return total mass source term [kg/m3/s]
                virtual tmp<fvScalarMatrix> Srho(volScalarField& rho) const = 0;


        // I-O

            //- Print cloud information
            virtual void info() = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
