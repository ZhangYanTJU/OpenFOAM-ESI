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

#include "SpaceChargeDensity.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::SpaceChargeDensity<CloudType>::write()
{
    if (erhoPtr_)
    {
        erhoPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "erhoPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SpaceChargeDensity<CloudType>::SpaceChargeDensity
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    mPtr_(nullptr),
    qPtr_(nullptr),
    cPtr_(nullptr),
    erhoPtr_(nullptr)
{}


template<class CloudType>
Foam::SpaceChargeDensity<CloudType>::SpaceChargeDensity
(
    const SpaceChargeDensity<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    mPtr_(nullptr),
    qPtr_(nullptr),
    cPtr_(nullptr),
    erhoPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SpaceChargeDensity<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    const fvMesh& mesh = this->owner().mesh();

    if (mPtr_)
    {
        mPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        mPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "m",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar(dimMass, Zero)
            )
        );
    }

    if (qPtr_)
    {
        qPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        qPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "q",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar(dimCurrent*dimTime, Zero)
            )
        );
    }

    if (cPtr_)
    {
        cPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        cPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "c",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar(dimMass*dimTime, Zero)
            )
        );
    }

    if (erhoPtr_)
    {
        erhoPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        erhoPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "erho",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimCurrent*dimTime/dimVolume, Zero)
            )
        );
    }
}


template<class CloudType>
void Foam::SpaceChargeDensity<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    volScalarField& m = mPtr_();
    volScalarField& q = qPtr_();
    volScalarField& c = cPtr_();
    volScalarField& erho = erhoPtr_();

    auto& own = this->owner();
    const fvMesh& mesh = own.mesh();

    c.primitiveFieldRef() /= mesh.time().deltaTValue()*mesh.V();

    // (YSSD:Eq. 4)

    // sum masses and charges in cells
    forAllConstIters(own, parcelIter)
    {
        const parcelType& p = parcelIter();
        m[p.cell()] += p.mass()*p.nParticle();
        q[p.cell()] += p.eq()*p.nParticle();
    }

    forAllConstIters(own, parcelIter)
    {
        const parcelType& p = parcelIter();

        // (YSSD:Eq. 5)
        erho[p.cell()] = q[p.cell()]/m[p.cell()]*c[p.cell()];
    }

    CloudFunctionObject<CloudType>::postEvolve(td);
}


template<class CloudType>
void Foam::SpaceChargeDensity<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point&,
    bool&
)
{
    volScalarField& c = cPtr_();

    // (YSSD:Eq. 3)
    c[p.cell()] += p.mass()*p.nParticle()*dt;
}


// ************************************************************************* //
