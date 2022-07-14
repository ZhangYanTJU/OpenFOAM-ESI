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

#include "fullIncompressibleVars.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fullIncompressibleVars, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void fullIncompressibleVars::calculateStorageMetrics()
{
    storageMetrics_.setSize(p_().storageMetrics().size(), Zero);
    label i = 0;
    if ( !storageParams_.timing() )
    {
        this->calculateAndWrite(p_(), pOld_, incoVars_.pInst().name(), i);
        this->calculateAndWrite(phi_(), phiOld_, incoVars_.phiInst().name(), i);
        this->calculateAndWrite(U_(), UOld_, incoVars_.UInst().name(), i);
        incompressible::RASModelVariables& rasVars =
            incoVars_.RASModelVariables()();
        if (rasVars.hasTMVar1())
        {
            this->calculateAndWrite
                (RASModelVars_[0], TMVar1Old_, rasVars.TMVar1Inst().name(), i);
        }
        if (rasVars.hasTMVar2())
        {
            this->calculateAndWrite
                (RASModelVars_[1], TMVar2Old_, rasVars.TMVar2Inst().name(), i);
        }
        storageMetrics_ = Zero;
    }
    addStorageMetricsContribution(p_(), pOld_);
    addStorageMetricsContribution(phi_(), phiOld_);
    addStorageMetricsContribution(U_(), UOld_);
    incompressible::RASModelVariables& rasVars =
        incoVars_.RASModelVariables()();
    if (rasVars.hasTMVar1())
    {
        addStorageMetricsContribution(RASModelVars_[0], TMVar1Old_);
    }
    if (rasVars.hasTMVar2())
    {
        addStorageMetricsContribution(RASModelVars_[1], TMVar2Old_);
    }
    if (!storageParams_.timing())
    {
        write(0);
        writeLog("All");
    }
    names_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fullIncompressibleVars::fullIncompressibleVars
(
    incompressibleVars& vs,
    storageParameters& storageParams,
    label a
)
:
    compressedIncompressibleVars(vs, storageParams),
    a_(a),
    pOld_(),
    phiOld_(),
    UOld_(),
    TMVar1Old_(),
    TMVar2Old_()
{
    if (a_ < 1)
    {
        FatalErrorInFunction
            << "Total number 'a' of time-solutions to be stored must be >= 1."
            << nl
            << "Object constructed with (a = " << a_ << ")"
            << endl
            << exit(FatalError);
    }
    storeOldTimes(incoVars_.pInst(), pOld_, 0);
    storeOldTimes(incoVars_.phiInst(), phiOld_, 1);
    storeOldTimes(incoVars_.UInst(), UOld_, 2);
    incompressible::RASModelVariables& rasVars =
        incoVars_.RASModelVariables()();
    if (rasVars.hasTMVar1())
    {
        storeOldTimes(rasVars.TMVar1Inst(), TMVar1Old_, 3);
    }
    if (rasVars.hasTMVar2())
    {
        storeOldTimes(rasVars.TMVar2Inst(), TMVar2Old_, 4);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fullIncompressibleVars::compress()
{
    //scalar startTime = mesh_.time().elapsedCpuTime();
    p_().compress();
    forAll(pOld_, iPtr)
    {
        pOld_[iPtr].compress();
    }
    phi_().compress();
    forAll(phiOld_, iPtr)
    {
        phiOld_[iPtr].compress();
    }
    U_().compress();
    forAll(UOld_, iPtr)
    {
        UOld_[iPtr].compress();
    }
    incompressible::RASModelVariables& rasVars =
        incoVars_.RASModelVariables()();
    if (rasVars.hasTMVar1())
    {
        RASModelVars_[0].compress();
        forAll(TMVar1Old_, iPtr)
        {
            TMVar1Old_[iPtr].compress();
        }
    }
    if (rasVars.hasTMVar2())
    {
        RASModelVars_[1].compress();
        forAll(TMVar2Old_, iPtr)
        {
            TMVar2Old_[iPtr].compress();
        }
    }
    //scalar endTime = mesh_.time().elapsedCpuTime();
    calculateStorageMetrics();
}


void fullIncompressibleVars::decompress
(
    incompressibleVars& vars
)
{
    // Phi must be retrieved first, since the BCs of p and U might depend on
    // flux
    decompressAll(vars.phiInst(), phi_(), phiOld_);
    DebugInfo
        << "DEBUG: phi field retrieved" << endl;
    decompressAll(vars.pInst(), p_(), pOld_);
    DebugInfo
        << "DEBUG: p field retrieved" << endl;
    decompressAll(vars.UInst(), U_(), UOld_);
    DebugInfo
        << "DEBUG: U field retrieved" << endl;
    incompressible::RASModelVariables& rasVars = vars.RASModelVariables()();
    if (rasVars.hasTMVar1())
    {
        decompressAll(rasVars.TMVar1Inst(), RASModelVars_[0], TMVar1Old_);
        DebugInfo
            << "DEBUG: 1st turbulence model field retrieved" << endl;
    }
    if (rasVars.hasTMVar2())
    {
        decompressAll(rasVars.TMVar2Inst(), RASModelVars_[1], TMVar2Old_);
        DebugInfo
            << "DEBUG: 2nd turbulence model field retrieved" << endl;
    }
}


void fullIncompressibleVars::decompress()
{
    this->decompress(incoVars_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
