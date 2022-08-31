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
    Foam::BaseCloud

Description
    Base template for parcel clouds. Inserts the baseCloudInterface
    virtualisation layer into the class. Also defines default zero-return
    source methods to enable the less functional clouds to be used in more
    complex situations.

\*---------------------------------------------------------------------------*/

#ifndef Foam_BaseCloud_H
#define Foam_BaseCloud_H

#include "Cloud.H"
#include "baseCloudInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class SLGThermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ParticleType>
class BaseCloud
:
    public Cloud<ParticleType>,
    virtual public baseCloudInterface
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;


public:

    // Constructors

        BaseCloud
        (
            const word& cloudName,
            const dimensionedVector& g,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const bool readFields = true
        )
        :
            Cloud<ParticleType>(rho.mesh(), cloudName, false),
            mesh_(rho.mesh())
        {}


        BaseCloud
        (
            const word& cloudName,
            const dimensionedVector& g,
            const volScalarField& rho,
            const volVectorField& U,
            const SLGThermo& thermo,
            const bool readFields = true
        )
        :
            Cloud<ParticleType>(rho.mesh(), cloudName, false),
            mesh_(rho.mesh())
        {}

        //- Copy constructor with new name
        BaseCloud
        (
            BaseCloud<ParticleType>& c,
            const word& name
        )
        :
            Cloud<ParticleType>(c.mesh(), name, c),
            mesh_(c.mesh_)
        {}

        //- Copy constructor with new name - creates bare cloud
        BaseCloud
        (
            const fvMesh& mesh,
            const word& name,
            const BaseCloud<ParticleType>& c
        )
        :
            Cloud<ParticleType>(mesh, name, IDLList<ParticleType>()),
            mesh_(mesh)
        {}

    //- Destructor
    virtual ~BaseCloud() = default;


    // Member Functions

        // Access

            //- Return reference to the mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }


        // Sources

            // Momentum

                //- Momentum transfer [kg m/s]
                virtual tmp<volVectorField::Internal> UTrans() const
                {
                    return volVectorField::Internal::New
                    (
                        IOobject::scopedName(this->name(), ":UTrans"),
                        mesh_,
                        dimensionedVector(dimMass*dimVelocity, Zero)
                    );
                }

                //- Momentum transfer coefficient [kg]
                virtual tmp<volScalarField::Internal> UCoeff() const
                {
                    return volScalarField::Internal::New
                    (
                        IOobject::scopedName(this->name(), "UCoeff"),
                        mesh_,
                        dimensionedScalar(dimMass, Zero)
                    );
                }

                //- Return momentum source term [kg m/s^2]
                virtual tmp<fvVectorMatrix> SU(volVectorField& U) const
                {
                    return tmp<fvVectorMatrix>::New(U, dimMass*dimAcceleration);
                }


            // Energy

                //- Return sensible enthalpy transfer [J]
                virtual tmp<volScalarField::Internal> hsTrans() const
                {
                    return volScalarField::Internal::New
                    (
                        IOobject::scopedName(this->name(), "hsTrans"),
                        mesh_,
                        dimensionedScalar(dimEnergy, Zero)
                    );
                }

                //- Return sensible enthalpy transfer coefficient [J/K]
                virtual tmp<volScalarField::Internal> hsCoeff() const
                {
                    return volScalarField::Internal::New
                    (
                        IOobject::scopedName(this->name(), "hsCoeff"),
                        mesh_,
                        dimensionedScalar(dimEnergy/dimTemperature, Zero)
                    );
                }

                //- Return sensible enthalpy source term [J/s]
                virtual tmp<fvScalarMatrix> Sh(volScalarField& hs) const
                {
                    return tmp<fvScalarMatrix>::New(hs, dimEnergy/dimTime);
                }


            // Radiation - overrides thermoCloud virtual abstract members

                //- Return equivalent particulate emission [kg/m/s^3]
                virtual tmp<volScalarField> Ep() const
                {
                    return volScalarField::New
                    (
                        IOobject::scopedName(this->name(), "radiation_Ep"),
                        mesh_,
                        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
                    );
                }

                //- Return equivalent particulate absorption [1/m]
                virtual tmp<volScalarField> ap() const
                {
                    return volScalarField::New
                    (
                        IOobject::scopedName(this->name(), "radiation_ap"),
                        mesh_,
                        dimensionedScalar(dimless/dimLength, Zero)
                    );
                }

                //- Return equivalent particulate scattering factor [1/m]
                virtual tmp<volScalarField> sigmap() const
                {
                    return volScalarField::New
                    (
                        IOobject::scopedName(this->name(), "radiation_sigmap"),
                        mesh_,
                        dimensionedScalar(dimless/dimLength, Zero)
                    );
                }


            // Mass

                //- Return mass source term for specie [kg/s]
                virtual tmp<fvScalarMatrix> SYi
                (
                    const label speciei,
                    volScalarField& Yi
                ) const
                {
                    return tmp<fvScalarMatrix>::New(Yi, dimMass/dimTime);
                }

                //- Return total mass source term [kg/m^3/s]
                virtual tmp<volScalarField::Internal> Srho() const
                {
                    return volScalarField::Internal::New
                    (
                        IOobject::scopedName(this->name(), "Srho"),
                        mesh_,
                        dimensionedScalar(dimDensity/dimTime, Zero)
                    );
                }

                //- Return total mass source term [kg/s]
                virtual tmp<fvScalarMatrix> Srho(volScalarField& rho) const
                {
                    return tmp<fvScalarMatrix>::New(rho, dimMass/dimTime);
                }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
