/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "meanVelocityForce.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(meanVelocityForce, 0);
    addToRunTimeSelectionTable(option, meanVelocityForce, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::meanVelocityForce::writeProps
(
    const scalar gradP
) const
{
    // Only write on output time
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::meanVelocityForce::meanVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    Ubar_(coeffs_.get<vector>("Ubar")),
    gradP0_(0.0),
    dGradP_(0.0),
    flowDir_(Ubar_/mag(Ubar_)),
    relaxation_(coeffs_.getOrDefault<scalar>("relaxation", 1)),
    rAPtr_(nullptr)
{
    coeffs_.readEntry("fields", fieldNames_);

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    fv::option::resetApplied();

    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timePath()/"uniform"/(name_ + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(propsFile);
        propsDict.readEntry("gradient", gradP0_);
    }

    Info<< "    Initial pressure gradient = " << gradP0_ << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::meanVelocityForce::magUbarAve
(
    const volVectorField& U
) const
{
    scalar magUbarAve = 0.0;

    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label celli = cells_[i];
        scalar volCell = cv[celli];
        magUbarAve += (flowDir_ & U[celli])*volCell;
    }

    reduce(magUbarAve, sumOp<scalar>());

    magUbarAve /= V_;

    return magUbarAve;
}


void Foam::fv::meanVelocityForce::correct(volVectorField& U)
{
    const scalarField& rAU = rAPtr_();

    // Integrate flow variables over cell set
    scalar rAUave = 0.0;
    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label celli = cells_[i];
        scalar volCell = cv[celli];
        rAUave += rAU[celli]*volCell;
    }

    // Collect across all processors
    reduce(rAUave, sumOp<scalar>());

    // Volume averages
    rAUave /= V_;

    scalar magUbarAve = this->magUbarAve(U);

    // Calculate the pressure gradient increment needed to adjust the average
    // flow-rate to the desired value
    dGradP_ = relaxation_*(mag(Ubar_) - magUbarAve)/rAUave;

    // Apply correction to velocity field
    forAll(cells_, i)
    {
        label celli = cells_[i];
        U[celli] += flowDir_*rAU[celli]*dGradP_;
    }

    U.correctBoundaryConditions();

    scalar gradP = gradP0_ + dGradP_;

    Info<< "Pressure gradient source: uncorrected Ubar = " << magUbarAve
        << ", pressure gradient = " << gradP << endl;

    writeProps(gradP);
}


void Foam::fv::meanVelocityForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    auto tSu = volVectorField::Internal::New
    (
        name_ + fieldNames_[fieldi] + "Sup",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );
    auto& Su = tSu.ref();

    scalar gradP = gradP0_ + dGradP_;

    UIndirectList<vector>(Su, cells_) = flowDir_*gradP;

    eqn += tSu;
}


void Foam::fv::meanVelocityForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    this->addSup(eqn, fieldi);
}


void Foam::fv::meanVelocityForce::constrain
(
    fvMatrix<vector>& eqn,
    const label
)
{
    if (!rAPtr_)
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                mesh_.newIOobject(IOobject::scopedName(name_, "rA")),
                1.0/eqn.A()
            )
        );
    }
    else
    {
        rAPtr_() = 1.0/eqn.A();
    }

    gradP0_ += dGradP_;
    dGradP_ = 0.0;
}


bool Foam::fv::meanVelocityForce::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}


// ************************************************************************* //
