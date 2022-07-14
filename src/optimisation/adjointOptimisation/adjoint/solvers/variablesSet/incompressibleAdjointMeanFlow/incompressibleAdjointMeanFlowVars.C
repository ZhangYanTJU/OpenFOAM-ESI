/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "incompressibleAdjointMeanFlowVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(incompressibleAdjointMeanFlowVars, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void incompressibleAdjointMeanFlowVars::setFields()
{
    setField(paPtr_, mesh_, "pa", solverName_, useSolverNameForFields_);
    setField(UaPtr_, mesh_, "Ua", solverName_, useSolverNameForFields_);
    setFluxField
    (
        phiaPtr_,
        mesh_,
        UaInst(),
        "phia",
        solverName_,
        useSolverNameForFields_
    );

    mesh_.setFluxRequired(paPtr_->name());
}


void incompressibleAdjointMeanFlowVars::setMeanFields()
{
    // Allocate mean fields
    // Only mean flow here since turbulent quantities
    // are allocated automatically in RASModelVariables
    if (solverControl_.average())
    {
        Info<< "Allocating Mean Adjoint Fields" << endl;
        variablesSet::setMeanField(paMeanPtr_, paInst(), mesh_);
        variablesSet::setMeanField(UaMeanPtr_, UaInst(), mesh_);
        variablesSet::setMeanField(phiaMeanPtr_, phiaInst(), mesh_);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

incompressibleAdjointMeanFlowVars::incompressibleAdjointMeanFlowVars
(
    fvMesh& mesh,
    solverControl& SolverControl,
    incompressibleVars& primalVars
)
:
    variablesSet(mesh, SolverControl.solverDict()),
    solverControl_(SolverControl),
    primalVars_(primalVars),
    paPtr_(nullptr),
    UaPtr_(nullptr),
    phiaPtr_(nullptr),
    paMeanPtr_(nullptr),
    UaMeanPtr_(nullptr),
    phiaMeanPtr_(nullptr)
{
    setFields();
    setMeanFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const incompressibleVars& incompressibleAdjointMeanFlowVars::primalVars() const
{
    return primalVars_;
}

const volScalarField& incompressibleAdjointMeanFlowVars::pa() const
{
    if (solverControl_.useAveragedFields())
    {
        return paMeanPtr_();
    }
    else
    {
        return paPtr_();
    }
}


volScalarField& incompressibleAdjointMeanFlowVars::pa()
{
    if (solverControl_.useAveragedFields())
    {
        return paMeanPtr_();
    }
    else
    {
        return paPtr_();
    }
}


const volVectorField& incompressibleAdjointMeanFlowVars::Ua() const
{
    if (solverControl_.useAveragedFields())
    {
        return UaMeanPtr_();
    }
    else
    {
        return UaPtr_();
    }
}


volVectorField& incompressibleAdjointMeanFlowVars::Ua()
{
    if (solverControl_.useAveragedFields())
    {
        return UaMeanPtr_();
    }
    else
    {
        return UaPtr_();
    }
}


const surfaceScalarField& incompressibleAdjointMeanFlowVars::phia() const
{
    if (solverControl_.useAveragedFields())
    {
        return phiaMeanPtr_();
    }
    else
    {
        return phiaPtr_();
    }
}


surfaceScalarField& incompressibleAdjointMeanFlowVars::phia()
{
    if (solverControl_.useAveragedFields())
    {
        return phiaMeanPtr_();
    }
    else
    {
        return phiaPtr_();
    }
}


bool incompressibleAdjointMeanFlowVars::computeMeanFields() const
{
    return solverControl_.average();
}


const solverControl& incompressibleAdjointMeanFlowVars::getSolverControl() const
{
    return solverControl_;
}


void incompressibleAdjointMeanFlowVars::nullify()
{
    // Nullify instantaneous fields
    variablesSet::nullifyField(paPtr_());
    variablesSet::nullifyField(UaPtr_());
    variablesSet::nullifyField(phiaPtr_());
}


void incompressibleAdjointMeanFlowVars::write()
{
    //  Write instantaneous fields
    variablesSet::writeField(paPtr_());
    variablesSet::writeField(UaPtr_());
    variablesSet::writeField(phiaPtr_());
    //  Write mean fields
    if (solverControl_.average())
    {
        variablesSet::writeField(paMeanPtr_());
        variablesSet::writeField(UaMeanPtr_());
        variablesSet::writeField(phiaMeanPtr_());
    }
}


void incompressibleAdjointMeanFlowVars::setWriteOption(IOobject::writeOption w)
{
    paPtr_().writeOpt() = w;
    UaPtr_().writeOpt() = w;
    phiaPtr_().writeOpt() = w;
    if (solverControl_.doAverageTime())
    {
        paMeanPtr_->writeOpt() = w;
        UaMeanPtr_->writeOpt() = w;
        phiaMeanPtr_->writeOpt() = w;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
