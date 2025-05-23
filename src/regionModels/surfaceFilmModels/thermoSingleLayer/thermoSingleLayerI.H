/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "thermoSingleLayer.H"
#include "filmRadiationModel.H"
#include "heatTransferModel.H"
#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const SLGThermo& thermoSingleLayer::thermo() const
{
    return thermo_;
}


inline tmp<scalarField> thermoSingleLayer::hs
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField& Cp = Cp_.boundaryField()[patchi];
    return Cp*(T - Tref.value());
}


inline tmp<volScalarField> thermoSingleLayer::hs
(
    const volScalarField& T
) const
{
    return volScalarField::New
    (
        "hs(" + T.name() + ")",
        IOobject::NO_REGISTER,
        Cp_*(T - Tref)
    );
}


inline tmp<volScalarField> thermoSingleLayer::T
(
    const volScalarField& hs
) const
{
    auto tT = volScalarField::New
    (
        "T(" + hs.name() + ")",
        IOobject::NO_REGISTER,
        hs/Cp_ + Tref
    );

    if (limitType::CLAMP_MIN == withTbounds_)
    {
        tT.ref().clamp_min(Tbounds_.min());
    }
    else if (limitType::CLAMP_MAX == withTbounds_)
    {
        tT.ref().clamp_max(Tbounds_.max());
    }
    else if (limitType::CLAMP_RANGE == withTbounds_)
    {
        tT.ref().clamp_range(Tbounds_);
    }

    return tT;
}


inline const volScalarField& thermoSingleLayer::hsSp() const
{
    return hsSp_;
}


inline const volScalarField& thermoSingleLayer::hsSpPrimary() const
{
    return hsSpPrimary_;
}


inline const volScalarField& thermoSingleLayer::TPrimary() const
{
    return TPrimary_;
}


inline const PtrList<volScalarField>& thermoSingleLayer::YPrimary() const
{
    return YPrimary_;
}


inline const heatTransferModel& thermoSingleLayer::htcs() const
{
    return *htcs_;
}


inline const heatTransferModel& thermoSingleLayer::htcw() const
{
    return *htcw_;
}


inline const phaseChangeModel& thermoSingleLayer::phaseChange() const
{
    return *phaseChange_;
}


inline const filmRadiationModel& thermoSingleLayer::radiation() const
{
    return *radiation_;
}


inline tmp<scalarField> thermoSingleLayer::qconvw(const label patchi) const
{
    const scalarField htc(htcw_->h()().boundaryField()[patchi]);
    const scalarField& Tp = T_.boundaryField()[patchi];
    const scalarField& Twp = Tw_.boundaryField()[patchi];

    return htc*(Tp - Twp);
}


inline tmp<scalarField> thermoSingleLayer::qconvp(const label patchi) const
{
    const scalarField htc(htcs_->h()().boundaryField()[patchi]);
    const scalarField& Tp = T_.boundaryField()[patchi];
    const scalarField& Tpp = TPrimary_.boundaryField()[patchi];

    return htc*(Tp - Tpp);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
