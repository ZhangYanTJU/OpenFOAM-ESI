/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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

#include "dynamicTopODesignVariables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicTopODesignVariables, 0);
    addToRunTimeSelectionTable
    (
        designVariables,
        dynamicTopODesignVariables,
        designVariables
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::dynamicTopODesignVariables::setActiveDesignVariables
(
    const label fluidID,
    const bool activeIO
)
{
    cellZoneMesh& cellZones = mesh_.cellZones();
    marchCells_.addFixedCells(cellZones, zones_.fixedPorousZoneIDs());
    marchCells_.addFixedCells(cellZones, zones_.fixedZeroPorousZoneIDs());
    if (!activeIO)
    {
        marchCells_.addFixedCells(zones_.IOCells());
    }

    marchCells_.update();
    activeDesignVariables_ = marchCells_.getActiveCells();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicTopODesignVariables::dynamicTopODesignVariables
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    dynamicTopODesignVariables(mesh, dict, mesh.nCells())
{}


Foam::dynamicTopODesignVariables::dynamicTopODesignVariables
(
    fvMesh& mesh,
    const dictionary& dict,
    const label size
)
:
    topODesignVariables(mesh, dict, size),
    marchCells_(mesh, dict.subDict("marchingCoeffs"))
{
    // Rest of the constructor initialization
    initialize();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::dynamicTopODesignVariables>
Foam::dynamicTopODesignVariables::New
(
    fvMesh& mesh,
    const dictionary& dict
)
{
    return autoPtr<dynamicTopODesignVariables>::New(mesh, dict);
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::dynamicTopODesignVariables::evolveNumber()
{
    marchCells_.update();
    activeDesignVariables_ = marchCells_.getActiveCells();
}


// ************************************************************************* //
