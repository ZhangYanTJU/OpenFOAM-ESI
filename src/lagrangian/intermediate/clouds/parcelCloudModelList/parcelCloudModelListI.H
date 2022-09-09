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

#include "fvMatrices.H"
#include "volFields.H"
#include "DimensionedField.H"

Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::parcelCloudModelList::UTrans() const
{
    auto tfld = volVectorField::Internal::New
    (
        "UTransEff",
        mesh_,
        dimensionedVector(dimMass*dimVelocity, Zero)
    );

    auto& fld = tfld.ref();

    for (const auto& c : *this)
    {
        fld += c.UTrans();
    }

    return tfld;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::parcelCloudModelList::UCeoff() const
{
    auto tfld = volScalarField::Internal::New
    (
        "UCoeffEff",
        mesh_,
        dimensionedScalar(dimMass, Zero)
    );

    auto& fld = tfld.ref();

    for (const auto& c : *this)
    {
        fld += c.UCoeff();
    }

    return tfld;
}


Foam::tmp<Foam::fvVectorMatrix> Foam::parcelCloudModelList::SU
(
    volVectorField& U
) const
{
    auto tfvm  = tmp<fvVectorMatrix>::New(U, dimForce);
    fvVectorMatrix& fvm = tfvm.ref();

    for (const auto& c : *this)
    {
        fvm += c.SU(U);
    }

    return tfvm;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::parcelCloudModelList::hsTrans() const
{
    auto tfld = volScalarField::Internal::New
    (
        "hsTransEff",
        mesh_,
        dimensionedScalar(dimEnergy, Zero)
    );

    auto& fld = tfld.ref();

    for (const auto& c : *this)
    {
        fld += c.hsTrans();
    }

    return tfld;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::parcelCloudModelList::hsCoeff() const
{
    auto tfld = volScalarField::Internal::New
    (
        "hsCoeffEff",
        mesh_,
        dimensionedScalar(dimEnergy/dimTemperature, Zero)
    );

    auto& fld = tfld.ref();

    for (const auto& c : *this)
    {
        fld += c.hsCoeff();
    }

    return tfld;
}

Foam::tmp<Foam::fvScalarMatrix> Foam::parcelCloudModelList::Sh
(
    volScalarField& hs
) const
{
    auto tfvm = tmp<fvScalarMatrix>::New(hs, dimEnergy/dimTime);
    fvScalarMatrix& fvm = tfvm.ref();

    for (const auto& c : *this)
    {
        fvm += c.Sh(hs);
    }

    return tfvm;
}


Foam::tmp<Foam::volScalarField> Foam::parcelCloudModelList::Ep() const
{
    auto tfld = volScalarField::New
    (
        "EpEff",
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    );

    auto& fld = tfld.ref();

    for (const auto& c : *this)
    {
        fld += c.Ep();
    }

    return tfld;
}


Foam::tmp<Foam::volScalarField> Foam::parcelCloudModelList::ap() const
{
    auto tfld = volScalarField::New
    (
        "apEff",
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    );

    auto& fld = tfld.ref();

    for (const auto& c : *this)
    {
        fld += c.ap();
    }

    return tfld;
}


Foam::tmp<Foam::volScalarField> Foam::parcelCloudModelList::sigmap() const
{
    auto tfld = volScalarField::New
    (
        "sigmapEff",
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    );

    auto& fld = tfld.ref();

    for (const auto& c : *this)
    {
        fld += c.sigmap();
    }

    return tfld;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::parcelCloudModelList::SYi
(
    const label speciei,
    volScalarField& Yi
) const
{
    auto tfvm = tmp<fvScalarMatrix>::New(Yi, dimMass/dimTime);
    fvScalarMatrix& fvm = tfvm.ref();

    for (const auto& c : *this)
    {
        fvm += c.SYi(speciei, Yi);
    }

    return tfvm;
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::parcelCloudModelList::Srho() const
{
    auto tfld = volScalarField::Internal::New
    (
        "SrhoEff",
        mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    );

    auto& fld = tfld.ref();

    for (const auto& c : *this)
    {
        fld += c.Srho();
    }

    return tfld;
}


Foam::tmp<Foam::fvScalarMatrix> Foam::parcelCloudModelList::Srho
(
    volScalarField& rho
) const
{
    auto tfvm = tmp<fvScalarMatrix>::New(rho, dimMass/dimTime);
    fvScalarMatrix& fvm = tfvm.ref();

    for (const auto& c : *this)
    {
        fvm += c.Srho(rho);
    }

    return tfvm;
}


// ************************************************************************* //
