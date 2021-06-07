/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "binModel.H"
#include "fvMesh.H"
#include "fluidThermo.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(binModel, 0);
    defineRunTimeSelectionTable(binModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binModel::binModel
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    porosity_(false),
    nBin_(1),
    patchSet_(),
    forceBinFilePtr_(nullptr),
    momentBinFilePtr_(nullptr),
    coeffBinFilePtrs_()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::binModel::read(const dictionary& dict)
{
    dict.readIfPresent("porosity", porosity_);

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        );

    if (!(dict.found("binModel")) && dict.found("binData")) // old behaviour
    {
        Info<< "    Selecting bin model: singleDirectionBin" << endl;
    }
    else
    {
        const word modelType(dict.getOrDefault<word>("binModel", "noBin"));
        Info<< "    Selecting bin model: " << modelType << endl;
    }

    coeffBinFilePtrs_.setSize(6);

    return true;
}


// ************************************************************************* //
