/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) 2020,2023 PCOpt/NTUA
    Copyright (C) 2020,2023 FOSS GP
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

#include "fvOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(option, 0);
        defineRunTimeSelectionTable(option, dictionary);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::option::resetApplied()
{
    applied_.resize_nocopy(fieldNames_.size());
    applied_ = false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::option::option
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    name_(name),
    modelType_(modelType),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.optionalSubDict(modelType + "Coeffs")),
    fieldNames_(),
    applied_(),
    active_(dict_.getOrDefault("active", true)),
    log(true)
{
    Log << incrIndent
        << indent << "Source: " << name_ << endl
        << indent << "State: " << (active_ ? "active" : "inactive") << endl
        << decrIndent;
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fv::option> Foam::fv::option::New
(
    const word& name,
    const dictionary& coeffs,
    const fvMesh& mesh
)
{
    const word modelType(coeffs.get<word>("type"));

    Info<< indent
        << "Selecting finite volume options type " << modelType << endl;

    mesh.time().libs().open
    (
        coeffs,
        "libs",
        dictionaryConstructorTablePtr_
    );

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            coeffs,
            "fvOption",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<fv::option>(ctorPtr(name, modelType, coeffs, mesh));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::option::isActive()
{
    return active_;
}


Foam::label Foam::fv::option::applyToField(const word& fieldName) const
{
    return fieldNames_.find(fieldName);
}


void Foam::fv::option::checkApplied() const
{
    forAll(applied_, i)
    {
        if (!applied_[i])
        {
            WarningInFunction
                << "Source " << name_ << " defined for field "
                << fieldNames_[i] << " but never used" << endl;
        }
    }
}


void Foam::fv::option::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<tensor>& eqn,
    const label fieldi
)
{
    addSup(alpha*rho, eqn, fieldi);
}


void Foam::fv::option::constrain(fvMatrix<scalar>& eqn, const label fieldi)
{}


void Foam::fv::option::constrain(fvMatrix<vector>& eqn, const label fieldi)
{}


void Foam::fv::option::constrain
(
    fvMatrix<sphericalTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::constrain
(
    fvMatrix<symmTensor>& eqn,
    const label fieldi
)
{}


void Foam::fv::option::constrain(fvMatrix<tensor>& eqn, const label fieldi)
{}


void Foam::fv::option::correct(volScalarField& field)
{}


void Foam::fv::option::correct(volVectorField& field)
{}


void Foam::fv::option::correct(volSphericalTensorField& field)
{}


void Foam::fv::option::correct(volSymmTensorField& field)
{}


void Foam::fv::option::correct(volTensorField& field)
{}


void Foam::fv::option::correct(surfaceScalarField& field)
{}


void Foam::fv::option::correct(surfaceVectorField& field)
{}


void Foam::fv::option::correct(surfaceTensorField& field)
{}


void Foam::fv::option::postProcessSens
(
    scalarField& sensField,
    const word& fieldName,
    const word& designVariablesName
)
{}


void Foam::fv::option::postProcessSens
(
    vectorField& sensField,
    const word& fieldName,
    const word& designVariablesName
)
{}


void Foam::fv::option::postProcessSens
(
    tensorField& sensField,
    const word& fieldName,
    const word& designVariablesName
)
{}


void Foam::fv::option::postProcessAuxSens
(
    const volScalarField& primalField,
    const volScalarField& adjointField,
    scalarField& sensField,
    const word& fieldName
)
{}


void Foam::fv::option::postProcessAuxSens
(
    const volVectorField& primalField,
    const volVectorField& adjointField,
    scalarField& sensField,
    const word& fieldName
)
{}


void Foam::fv::option::postProcessAuxSens
(
    const volTensorField& primalField,
    const volTensorField& adjointField,
    scalarField& sensField,
    const word& fieldName
)
{}


// ************************************************************************* //
