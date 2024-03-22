/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "treeTurbulence.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(treeTurbulence, 0);
    addToRunTimeSelectionTable(option, treeTurbulence, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::volScalarField& Foam::fv::treeTurbulence::getOrReadField
(
    const word& fieldName
) const
{
    auto* ptr = mesh_.getObjectPtr<volScalarField>(fieldName);

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        );
        mesh_.objectRegistry::store(ptr);
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::treeTurbulence::treeTurbulence
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    isEpsilon_(true),
    betaP_(),
    betaD_(),
    Ceps1_(),
    Ceps2_(),
    betaStar_(),
    CdName_(),
    LADname_()
{
    read(dict);

    const auto* turbPtr = mesh_.findObject<turbulenceModel>
    (
        turbulenceModel::propertiesName
    );

    if (!turbPtr)
    {
        FatalErrorInFunction
            << "Unable to find a turbulence model."
            << abort(FatalError);
    }

    fieldNames_.resize(2);

    tmp<volScalarField> tepsilon = turbPtr->epsilon();
    tmp<volScalarField> tomega = turbPtr->omega();

    if (!tepsilon.isTmp())
    {
        fieldNames_[0] = tepsilon().name();
    }
    else if (!tomega.isTmp())
    {
        isEpsilon_ = false;
        fieldNames_[0] = tomega().name();
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find epsilon or omega field." << nl
            << abort(FatalError);
    }

    fieldNames_[1] = turbPtr->k()().name();

    fv::option::resetApplied();

    Log << "    Applying treeTurbulence to: "
        << fieldNames_[0] << " and " << fieldNames_[1]
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::treeTurbulence::addSup
(
    fvScalarMatrix& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        kSource(geometricOneField(), geometricOneField(), eqn, fieldi);
        return;
    }

    if (isEpsilon_)
    {
        epsilonSource(geometricOneField(), geometricOneField(), eqn, fieldi);
    }
    else
    {
        omegaSource(geometricOneField(), geometricOneField(), eqn, fieldi);
    }
}


void Foam::fv::treeTurbulence::addSup
(
    const volScalarField& rho,
    fvScalarMatrix& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        kSource(geometricOneField(), rho, eqn, fieldi);
        return;
    }

    if (isEpsilon_)
    {
        epsilonSource(geometricOneField(), rho, eqn, fieldi);
    }
    else
    {
        omegaSource(geometricOneField(), rho, eqn, fieldi);
    }
}


void Foam::fv::treeTurbulence::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvScalarMatrix& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        kSource(alpha, rho, eqn, fieldi);
        return;
    }

    if (isEpsilon_)
    {
        epsilonSource(alpha, rho, eqn, fieldi);
    }
    else
    {
        omegaSource(alpha, rho, eqn, fieldi);
    }
}


bool Foam::fv::treeTurbulence::read(const dictionary& dict)
{
    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }

    betaP_ = dict.getOrDefault<scalar>("betaP", 1.0);
    betaD_ = dict.getOrDefault<scalar>("betaD", 4.0);
    Ceps1_ = dict.getOrDefault<scalar>("Ceps1", 0.9);
    Ceps2_ = dict.getOrDefault<scalar>("Ceps2", 0.9);
    betaStar_ = dict.getOrDefault<scalar>("betaStar", 0.09);
    CdName_ = dict.getOrDefault<word>("Cd", "Cd");
    LADname_ = dict.getOrDefault<word>("LAD", "LAD");

    (void) getOrReadField(CdName_);
    (void) getOrReadField(LADname_);

    return true;
}


// ************************************************************************* //
