/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "MultiComponentPhaseModel.H"

#include "phaseSystem.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MultiComponentPhaseModel<BasePhaseModel>::MultiComponentPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
:
    BasePhaseModel(fluid, phaseName, index),
    Sct_
    (
        "Sct",
        dimless,
        fluid.subDict(phaseName)
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        fluid.mesh().solverDict("Yi")
    ),
    inertIndex_(-1)
{
    word inertSpecie;
    if
    (
        this->thermo_->readIfPresent("inertSpecie", inertSpecie)
     && !inertSpecie.empty()
    )
    {
        inertIndex_ = this->thermo_->composition().species().find(inertSpecie);
    }

    PtrList<volScalarField>& Y = this->thermo_->composition().Y();

    forAll(Y, i)
    {
        if (i != inertIndex_ && this->thermo_->composition().active(i))
        {
            const label j = YActive_.size();
            YActive_.resize(j + 1);
            YActive_.set(j, &Y[i]);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MultiComponentPhaseModel<BasePhaseModel>::~MultiComponentPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MultiComponentPhaseModel<BasePhaseModel>::correctThermo()
{
    volScalarField Yt
    (
        IOobject
        (
            IOobject::groupName("Yt", this->name()),
            this->fluid().mesh().time().timeName(),
            this->fluid().mesh()
        ),
        this->fluid().mesh(),
        dimensionedScalar(dimless)
    );

    PtrList<volScalarField>& Yi = YRef();

    forAll(Yi, i)
    {
        if (i != inertIndex_)
        {
            Yt += Yi[i];
        }
    }

    if (inertIndex_ != -1)
    {
        Yi[inertIndex_] = scalar(1) - Yt;
        Yi[inertIndex_].clamp_min(0);
    }
    else
    {
        forAll(Yi, i)
        {
            Yi[i] /= Yt;
            Yi[i].clamp_min(0);
        }
    }

    BasePhaseModel::correctThermo();
}


template<class BasePhaseModel>
bool Foam::MultiComponentPhaseModel<BasePhaseModel>::pure() const
{
    return false;
}


template<class BasePhaseModel>
Foam::tmp<Foam::fvScalarMatrix>
Foam::MultiComponentPhaseModel<BasePhaseModel>::YiEqn(volScalarField& Yi)
{
    const volScalarField& alpha = *this;
    const surfaceScalarField alphaRhoPhi(this->alphaRhoPhi());
    const tmp<volScalarField> trho(this->thermo().rho());
    const volScalarField& rho = trho();

    return
    (
        fvm::ddt(alpha, rho, Yi)
      + fvm::div(alphaRhoPhi, Yi, "div(" + alphaRhoPhi.name() + ",Yi)")

      - fvm::laplacian
        (
            fvc::interpolate(alpha)
           *fvc::interpolate(this->muEff()/Sct_),
            Yi
        )
     ==
        alpha*this->R(Yi)

      + fvc::ddt(residualAlpha_*rho, Yi)
      - fvm::ddt(residualAlpha_*rho, Yi)
    );
}


template<class BasePhaseModel>
const Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::Y() const
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
const Foam::volScalarField&
Foam::MultiComponentPhaseModel<BasePhaseModel>::Y(const word& name) const
{
    return this->thermo_->composition().Y(name);
}


template<class BasePhaseModel>
Foam::PtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::YRef()
{
    return this->thermo_->composition().Y();
}


template<class BasePhaseModel>
const Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::YActive() const
{
    return YActive_;
}


template<class BasePhaseModel>
Foam::UPtrList<Foam::volScalarField>&
Foam::MultiComponentPhaseModel<BasePhaseModel>::YActiveRef()
{
    return YActive_;
}


// ************************************************************************* //
