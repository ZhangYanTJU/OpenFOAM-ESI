/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#include "thermalShellModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalShellModel, 0);

defineRunTimeSelectionTable(thermalShellModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalShellModel::thermalShellModel
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    regionFaModel(mesh, "thermalShell", modelType, dict, true),
    TName_(dict.getOrDefault<word>("T", suffixed("Ts"))),
    TprimaryName_(dict.getOrDefault<word>("Tprimary", "T")),
    Tp_(mesh.lookupObject<volScalarField>(TprimaryName_)),
    T_
    (
        IOobject
        (
            TName_,
            regionMesh().time().timeName(),
            regionMesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    faOptions_
    (
        Foam::fa::options::New(primaryMesh(), regionFaModel::areaName())
    )
{
    if (faOptions_.optionList::empty())
    {
        Info<< "No finite-area options present for area:"
            << regionFaModel::areaName() << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermalShellModel::preEvolveRegion()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
