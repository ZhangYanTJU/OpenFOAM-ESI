/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "ThermoCloud.H"
#include "ThermoParcel.H"

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoCloud<CloudType>::setModels()
{
    heatTransferModel_.reset
    (
        HeatTransferModel<ThermoCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    TIntegrator_.reset
    (
        integrationScheme::New
        (
            "T",
            this->solution().integrationSchemes()
        ).ptr()
    );

    this->subModelProperties().readEntry("radiation", radiation_);

    if (radiation_)
    {
        radAreaP_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    IOobject::scopedName(this->name(), "radAreaP"),
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                this->mesh(),
                dimensionedScalar(dimArea, Zero)
            )
        );

        radT4_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    IOobject::scopedName(this->name(), "radT4"),
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                this->mesh(),
                dimensionedScalar(pow4(dimTemperature), Zero)
            )
        );

        radAreaPT4_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    IOobject::scopedName(this->name(), "radAreaPT4"),
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                this->mesh(),
                dimensionedScalar(sqr(dimLength)*pow4(dimTemperature), Zero)
            )
        );
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::cloudReset(ThermoCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    heatTransferModel_.reset(c.heatTransferModel_.ptr());
    TIntegrator_.reset(c.TIntegrator_.ptr());

    radiation_ = c.radiation_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType
    (
        cloudName,
        rho,
        U,
        thermo.thermo().mu(),
        g,
        false
    ),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(this->particleProperties()),
    thermo_(thermo),
    T_(thermo.thermo().T()),
    p_(thermo.thermo().p()),
    heatTransferModel_(nullptr),
    TIntegrator_(nullptr),
    radiation_(false),
    radAreaP_(nullptr),
    radT4_(nullptr),
    radAreaPT4_(nullptr),
    hsTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::scopedName(this->name(), "hsTrans"),
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy, Zero)
        )
    ),
    hsCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::scopedName(this->name(), "hsCoeff"),
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            this->mesh(),
            dimensionedScalar(dimEnergy/dimTemperature, Zero)
        )
    )
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
            this->deleteLostParticles();
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    ThermoCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(c.constProps_),
    thermo_(c.thermo_),
    T_(c.T()),
    p_(c.p()),
    heatTransferModel_(c.heatTransferModel_->clone()),
    TIntegrator_(c.TIntegrator_->clone()),
    radiation_(c.radiation_),
    radAreaP_(nullptr),
    radT4_(nullptr),
    radAreaPT4_(nullptr),
    hsTrans_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::scopedName(this->name(), "hsTrans"),
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            c.hsTrans()
        )
    ),
    hsCoeff_
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::scopedName(this->name(), "hsCoeff"),
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            c.hsCoeff()
        )
    )
{
    if (radiation_)
    {
        radAreaP_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    IOobject::scopedName(this->name(), "radAreaP"),
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                c.radAreaP()
            )
        );

        radT4_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    IOobject::scopedName(this->name(), "radT4"),
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                c.radT4()
            )
        );

        radAreaPT4_.reset
        (
            new volScalarField::Internal
            (
                IOobject
                (
                    IOobject::scopedName(this->name(), "radAreaPT4"),
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                c.radAreaPT4()
            )
        );
    }
}


template<class CloudType>
Foam::ThermoCloud<CloudType>::ThermoCloud
(
    const fvMesh& mesh,
    const word& name,
    const ThermoCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    thermoCloud(),
    cloudCopyPtr_(nullptr),
    constProps_(),
    thermo_(c.thermo()),
    T_(c.T()),
    p_(c.p()),
    heatTransferModel_(nullptr),
    TIntegrator_(nullptr),
    radiation_(false),
    radAreaP_(nullptr),
    radT4_(nullptr),
    radAreaPT4_(nullptr),
    hsTrans_(nullptr),
    hsCoeff_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ThermoCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    parcel.T() = constProps_.T0();
    parcel.Cp() = constProps_.Cp0();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ThermoCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    hsTrans_->field() = 0.0;
    hsCoeff_->field() = 0.0;

    if (radiation_)
    {
        radAreaP_->field() = 0.0;
        radT4_->field() = 0.0;
        radAreaPT4_->field() = 0.0;
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::relaxSources
(
    const ThermoCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    this->relax(hsTrans_(), cloudOldTime.hsTrans(), "h");
    this->relax(hsCoeff_(), cloudOldTime.hsCoeff(), "h");

    if (radiation_)
    {
        this->relax(radAreaP_(), cloudOldTime.radAreaP(), "radiation");
        this->relax(radT4_(), cloudOldTime.radT4(), "radiation");
        this->relax(radAreaPT4_(), cloudOldTime.radAreaPT4(), "radiation");
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    this->scale(hsTrans_(), "h");
    this->scale(hsCoeff_(), "h");

    if (radiation_)
    {
        this->scale(radAreaP_(), "radiation");
        this->scale(radT4_(), "radiation");
        this->scale(radAreaPT4_(), "radiation");
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{
    CloudType::preEvolve(td);

    this->pAmbient() = thermo_.thermo().p().average().value();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ThermoCloud<CloudType>::info()
{
    CloudType::info();

    Log_<< "    Temperature min/max             = " << Tmin() << ", " << Tmax()
        << endl;
}


// ************************************************************************* //
