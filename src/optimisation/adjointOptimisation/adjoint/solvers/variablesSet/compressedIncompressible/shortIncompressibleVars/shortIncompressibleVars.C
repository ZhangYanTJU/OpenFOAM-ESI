/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 PCOpt/NTUA
    Copyright (C) 2022      FOSS GP
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

#include "shortIncompressibleVars.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"
#include "profiling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(shortIncompressibleVars, 0);
addToRunTimeSelectionTable
(
    compressedIncompressibleVars,
    shortIncompressibleVars,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

shortIncompressibleVars::shortIncompressibleVars
(
    incompressibleVars& vs,
    storageParameters& storageParams
)
:
    compressedIncompressibleVars(vs, storageParams)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void shortIncompressibleVars::compress()
{
    addProfiling(shortIncompressibleVars, "shortIncompressibleVars::compress");
    //scalar startTime = mesh_.time().elapsedCpuTime();
    p_().compress();
    DebugInfo
        << "DEBUG: p field stored" << endl;
    phi_().compress();
    DebugInfo
        << "DEBUG: phi field stored" << endl;
    U_().compress();
    DebugInfo
        << "DEBUG: U field stored" << endl;
    forAll (RASModelVars_, iPtr)
    {
        RASModelVars_[iPtr].compress();
        DebugInfo
            << "DEBUG: " << iPtr + 1 << "th turbulence model field stored"
            << endl;
    }
    //scalar endTime = mesh_.time().elapsedCpuTime();
}


void shortIncompressibleVars::decompress
(
    incompressibleVars& vars
)
{
    addProfiling
        (shortIncompressibleVars, "shortIncompressibleVars::decompress");
    // Phi must be retrieved first, since the BCs of p and U might depend on
    // flux
    phi_().decompress(vars.phiInst());
    DebugInfo
        << "DEBUG: phi field retrieved" << endl;
    p_().decompress(vars.pInst());
    DebugInfo
        << "DEBUG: p field retrieved" << endl;
    U_().decompress(vars.UInst());
    DebugInfo
        << "DEBUG: U field retrieved" << endl;
    incompressible::RASModelVariables& rasVars = vars.RASModelVariables()();
    if (rasVars.hasTMVar1())
    {
        RASModelVars_[0].decompress(rasVars.TMVar1Inst());
        DebugInfo
            << "DEBUG: 1st turbulence model field retrieved" << endl;
    }
    if (rasVars.hasTMVar2())
    {
        RASModelVars_[1].decompress(rasVars.TMVar2Inst());
        DebugInfo
            << "DEBUG: 2nd turbulence model field retrieved" << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
