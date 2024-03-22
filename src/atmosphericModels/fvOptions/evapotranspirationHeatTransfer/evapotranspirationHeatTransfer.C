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

#include "evapotranspirationHeatTransfer.H"
#include "evapotranspirationHeatTransferModel.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(evapotranspirationHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        option,
        evapotranspirationHeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::evapotranspirationHeatTransfer::evapotranspirationHeatTransfer
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    ethtModelPtr_()
{
    // Set the field name to that of the energy
    // field from which the temperature is obtained

    const auto& thermo = mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.resize(1, thermo.he().name());

    fv::option::resetApplied();

    Info<< "    Applying evapotranspirationHeatTransfer to: " << fieldNames_[0]
        << endl;

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::evapotranspirationHeatTransfer::~evapotranspirationHeatTransfer()
{}  // evapotranspirationHeatTransferModel was forward declared


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::evapotranspirationHeatTransfer::addSup
(
    fvScalarMatrix& eqn,
    const label fieldi
)
{
    if (this->V() < VSMALL)
    {
        return;
    }

    const scalar V = this->V();
    const scalarField& Vcells = mesh_.V();

    const volScalarField& he = eqn.psi();
    scalarField& heSource = eqn.source();

    // Calculate evapotranspiration heat transfer rate per volume [J/s/m^3]
    const scalarField Q(ethtModelPtr_->Q(cells_)/V);

    if (he.dimensions() == dimEnergy/dimMass)
    {
        forAll(cells_, i)
        {
            const label celli = cells_[i];

            heSource[celli] += Q[i]*Vcells[celli];
        }
    }
    else if (he.dimensions() == dimTemperature)
    {
        const auto& thermo =
            mesh_.lookupObject<basicThermo>(basicThermo::dictName);

        // Calculate density*heat capacity at constant pressure/volume
        const volScalarField rhoCpv(thermo.rho()*thermo.Cpv());

        // heSource [K m^3/s] = [J/s/m^3] * m^3 / [kg/m^3] / [J/kg/K]
        forAll(cells_, i)
        {
            const label celli = cells_[i];

            heSource[celli] += Q[i]*Vcells[celli]*rhoCpv[i];
        }
    }
}


void Foam::fv::evapotranspirationHeatTransfer::addSup
(
    const volScalarField& rho,
    fvScalarMatrix& eqn,
    const label fieldi
)
{
    addSup(eqn, fieldi);
}


bool Foam::fv::evapotranspirationHeatTransfer::read(const dictionary& dict)
{
    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }

    if (selectionMode_ != selectionModeType::smCellZone)
    {
        FatalIOErrorInFunction(dict)
            << "evapotranspirationHeatTransfer requires "
            << selectionModeTypeNames_[selectionModeType::smCellZone]
            << exit(FatalIOError);
    }

    ethtModelPtr_.reset
    (
        evapotranspirationHeatTransferModel::New(dict, mesh_)
    );

    return true;
}


// ************************************************************************* //
