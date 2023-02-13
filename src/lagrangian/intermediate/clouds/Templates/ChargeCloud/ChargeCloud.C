/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "ChargeCloud.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ChargeCloud<CloudType>::setModels()
{
}


template<class CloudType>
void Foam::ChargeCloud<CloudType>::cloudReset
(
    ChargeCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ChargeCloud<CloudType>::ChargeCloud
(
    const word& cloudName,
    const dimensionedVector& g,
    const volScalarField& rho,
    const volVectorField& U,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType(cloudName, g, rho, U, thermo, false),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties())
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
            this->deleteLostParticles();
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ChargeCloud<CloudType>::ChargeCloud
(
    ChargeCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    cloudCopyPtr_(nullptr),
    constProps_(c.constProps_)
{}


template<class CloudType>
Foam::ChargeCloud<CloudType>::ChargeCloud
(
    const fvMesh& mesh,
    const word& name,
    const ChargeCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    cloudCopyPtr_(nullptr),
    constProps_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ChargeCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);
}


template<class CloudType>
void Foam::ChargeCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    if (fullyDescribed)
    {

    }
}


template<class CloudType>
void Foam::ChargeCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ChargeCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ChargeCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ChargeCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
}


template<class CloudType>
void Foam::ChargeCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::ChargeCloud<CloudType>::autoMap
(
    const mapPolyMesh& mapper
)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ChargeCloud<CloudType>::info()
{
    CloudType::info();
}

/*
template<class CloudType>
void Foam::ChargeCloud<CloudType>::writeFields() const
{
    if (this->compositionModel_)
    {
        CloudType::particleType::writeFields(*this, this->composition());
    }
}
*/

// ************************************************************************* //
