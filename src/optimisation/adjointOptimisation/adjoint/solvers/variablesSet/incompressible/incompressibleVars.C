/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "incompressibleVars.H"
#include "createZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(incompressibleVars, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void incompressibleVars::setFields()
{
    setField(pPtr_, mesh_, "p", solverName_, useSolverNameForFields_);
    setField(UPtr_, mesh_, "U", solverName_, useSolverNameForFields_);
    setFluxField
    (
        phiPtr_,
        mesh_,
        UInst(),
        "phi",
        solverName_,
        useSolverNameForFields_
    );

    mesh_.setFluxRequired(pPtr_->name());

    // if required, correct boundary conditions of mean flow fields here in
    // order to have the correct bcs for e.g. turbulence models that follow.
    // NOTE: phi correction depends on the solver (includes for instance
    // Rhie-Chow interpolation).  This needs to be implemented within
    // incompressiblePrimalSolver
    if (correctBoundaryConditions_)
    {
        correctNonTurbulentBoundaryConditions();
    }

    laminarTransportPtr_.reset
    (
        new singlePhaseTransportModel(UInst(), phiInst())
    );
    turbulence_.reset
    (
        incompressible::turbulenceModel::New
        (
            UInst(),
            phiInst(),
            laminarTransport()
        ).ptr()
    );
    RASModelVariables_.reset
    (
        incompressible::RASModelVariables::New
        (
            mesh_,
            solverControl_
        ).ptr()
    );
    renameTurbulenceFields();
    if (correctBoundaryConditions_)
    {
        correctTurbulentBoundaryConditions();
    }
}


void incompressibleVars::setInitFields()
{
    // Store init fields
    // only mean flow here since turbulent quantities
    // are allocated automatically in RASModelVariables
    if (solverControl_.storeInitValues())
    {
        DebugInfo
            << "Allocating Initial Primal Fields" << endl;
        setInitField(pInitPtr_, pInst());
        setInitField(UInitPtr_, UInst());
        setInitField(phiInitPtr_, phiInst());
    }
}


void incompressibleVars::setMeanFields()
{
    // Allocate mean fields
    // only mean flow here since turbulent quantities
    // are allocated automatically in RASModelVariables
    if (solverControl_.average())
    {
        Info<< "Allocating Mean Primal Fields" << endl;
        setMeanField(pMeanPtr_, pInst(), mesh_);
        setMeanField(UMeanPtr_, UInst(), mesh_);
        setMeanField(phiMeanPtr_, phiInst(), mesh_);
        // Correct boundary conditions if necessary
        if (correctBoundaryConditions_)
        {
            pMeanPtr_().correctBoundaryConditions();
            UMeanPtr_().correctBoundaryConditions();
        }
    }
}


void incompressibleVars::renameTurbulenceFields()
{
    //  Turbulence model always reads fields with the prescribed name
    //  If a custom name is supplied, check whether this field exists,
    //  copy it to the field known by the turbulence model
    //  and re-name the latter
    if (useSolverNameForFields_)
    {
        incompressible::RASModelVariables& rasVars = RASModelVariables_();
        if (rasVars.hasTMVar1())
        {
            renameTurbulenceField(rasVars.TMVar1Inst(), solverName_);
        }
        if (rasVars.hasTMVar2())
        {
            renameTurbulenceField(rasVars.TMVar2Inst(), solverName_);
        }
        if (rasVars.hasNut())
        {
            renameTurbulenceField(rasVars.nutRefInst(), solverName_);
        }
    }
}


void incompressibleVars::correctNonTurbulentBoundaryConditions()
{
    Info<< "Correcting (U,p) boundary conditions " << endl;
    pInst().correctBoundaryConditions();
    UInst().correctBoundaryConditions();
    if (solverControl_.average())
    {
        pMeanPtr_().correctBoundaryConditions();
        UMeanPtr_().correctBoundaryConditions();
    }
}


void incompressibleVars::correctTurbulentBoundaryConditions()
{
    // If required, correct boundary conditions of turbulent fields.
    // Includes the correction of boundary conditions for averaged fields,
    // if any
    Info<< "Correcting boundary conditions of turbulent fields" << endl;
    RASModelVariables_().correctBoundaryConditions(turbulence_());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

incompressibleVars::incompressibleVars
(
    fvMesh& mesh,
    solverControl& SolverControl,
    scalar readTime
)
:
    variablesSet(mesh, SolverControl.solverDict()),
    solverControl_(SolverControl),
    pPtr_(nullptr),
    UPtr_(nullptr),
    phiPtr_(nullptr),
    laminarTransportPtr_(nullptr),
    turbulence_(nullptr),
    RASModelVariables_(nullptr),

    pInitPtr_(),
    UInitPtr_(),
    phiInitPtr_(),

    pMeanPtr_(nullptr),
    UMeanPtr_(nullptr),
    phiMeanPtr_(nullptr),

    correctBoundaryConditions_
    (
        SolverControl.solverDict().subOrEmptyDict("fieldReconstruction").
            getOrDefault<bool>("reconstruct", false)
    ),
    solDirs_((mesh.solutionD() + Vector<label>::one)/2),
    writeFields_
    (
        SolverControl.solverDict().getOrDefault<bool>("writeFields", true)
    )
{
    // Necessary for continuation of unsteady solvers. For steady solvers
    // should have no effect
    Time& time = const_cast<Time&>(mesh_.time());
    const dimensionedScalar startTime = mesh_.time();
    const label startTimeIndex = mesh_.time().timeIndex();
    if (readTime != -1)
    {
        // Set time to readTime. Time-index remains unchanged, since
        // it does not influence fields reading
        DebugInfo
            << "Initializing primal fields begin: time changes to "
            << ::Foam::name(readTime) << endl;
        time.setTime(readTime, startTimeIndex);
    }
    setFields();
    setInitFields();
    setMeanFields();
    storeAllocatedFieldNames();
    if (readTime != -1)
    {
        DebugInfo
            << "Initializing primal fields end: time changes to "
            << ::Foam::name(startTime.value()) << endl;
        // Restore time
        time.setTime(startTime, startTimeIndex);
    }
    if (!writeFields_)
    {
        setWriteOption(IOobject::NO_WRITE);
    }
}


incompressibleVars::incompressibleVars
(
    const incompressibleVars& vs
)
:
    variablesSet(vs.mesh_, vs.solverControl_.solverDict()),
    solverControl_(vs.solverControl_),
    pPtr_(allocateRenamedField(vs.pPtr_)),
    UPtr_(allocateRenamedField(vs.UPtr_)),
    phiPtr_(allocateRenamedField(vs.phiPtr_)),
    laminarTransportPtr_(nullptr),
    turbulence_(nullptr),
    RASModelVariables_(vs.RASModelVariables_.clone()),

    pInitPtr_(vs.pInitPtr_.size()),
    UInitPtr_(vs.UInitPtr_.size()),
    phiInitPtr_(vs.phiInitPtr_.size()),

    pMeanPtr_(allocateRenamedField(vs.pMeanPtr_)),
    UMeanPtr_(allocateRenamedField(vs.UMeanPtr_)),
    phiMeanPtr_(allocateRenamedField(vs.phiMeanPtr_)),

    correctBoundaryConditions_(vs.correctBoundaryConditions_),
    allocatedFieldNames_(vs.allocatedFieldNames_),
    solDirs_(vs.solDirs_)
{
    DebugInfo
        << "Calling incompressibleVars copy constructor" << endl;
    forAll(pInitPtr_, i)
    {
        pInitPtr_.set(i, allocateRenamedField(vs.pInitPtr_[i]));
    }
    forAll(UInitPtr_, i)
    {
        UInitPtr_.set(i, allocateRenamedField(vs.UInitPtr_[i]));
    }
    forAll(phiInitPtr_, i)
    {
        phiInitPtr_.set(i, allocateRenamedField(vs.phiInitPtr_[i]));
    }
}


autoPtr<variablesSet> incompressibleVars::clone() const
{
    DebugInfo
        << "Calling incompressibleVars::clone" << endl;

    return autoPtr<variablesSet>(new incompressibleVars(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void incompressibleVars::storeAllocatedFieldNames()
{
    allocatedFieldNames_.setSize(3);
    allocatedFieldNames_[0] = pInst().name();
    allocatedFieldNames_[1] = phiInst().name();
    allocatedFieldNames_[2] = UInst().name();
    if (RASModelVariables_().hasTMVar1())
    {
        allocatedFieldNames_.append(RASModelVariables_().TMVar1Inst().name());
    }
    if (RASModelVariables_().hasTMVar2())
    {
        allocatedFieldNames_.append(RASModelVariables_().TMVar2Inst().name());
    }
}


const volScalarField& incompressibleVars::p() const
{
    if (solverControl_.useAveragedFields())
    {
        return pMeanPtr_();
    }
    else
    {
        return pPtr_();
    }
}


volScalarField& incompressibleVars::p()
{
    if (solverControl_.useAveragedFields())
    {
        return pMeanPtr_();
    }
    else
    {
        return pPtr_();
    }
}


const volVectorField& incompressibleVars::U() const
{
    if (solverControl_.useAveragedFields())
    {
        return UMeanPtr_();
    }
    else
    {
        return UPtr_();
    }
}


volVectorField& incompressibleVars::U()
{
    if (solverControl_.useAveragedFields())
    {
        return UMeanPtr_();
    }
    else
    {
        return UPtr_();
    }
}


const surfaceScalarField& incompressibleVars::phi() const
{
    if (solverControl_.useAveragedFields())
    {
        return phiMeanPtr_();
    }
    else
    {
        return phiPtr_();
    }
}


surfaceScalarField& incompressibleVars::phi()
{
    if (solverControl_.useAveragedFields())
    {
        return phiMeanPtr_();
    }
    else
    {
        return phiPtr_();
    }
}


void incompressibleVars::restoreInitValues()
{
    if (solverControl_.storeInitValues())
    {
        Info<< "Restoring field values to initial ones" << endl;
        variablesSet::restoreFieldInitialization(pInitPtr_, pInst());
        variablesSet::restoreFieldInitialization(phiInitPtr_, phiInst());
        variablesSet::restoreFieldInitialization(UInitPtr_, UInst());
        RASModelVariables_().restoreInitValues();
    }
}


void incompressibleVars::resetMeanFields()
{
    if (solverControl_.average())
    {
        Info<< "Resetting mean fields to zero" << endl;

        // Reset fields to zero
        pMeanPtr_() == dimensionedScalar(pInst().dimensions(), Zero);
        UMeanPtr_() == dimensionedVector(UInst().dimensions(), Zero);
        phiMeanPtr_() == dimensionedScalar(phiInst().dimensions(), Zero);
        RASModelVariables_().resetMeanFields();

        // Reset averaging iteration index to 0
        solverControl_.averageIter() = 0;
    }
}


void incompressibleVars::computeMeanFields()
{
    if (solverControl_.doAverageIter())
    {
        Info<< "Averaging fields" << endl;
        label& iAverageIter = solverControl_.averageIter();
        scalar avIter(iAverageIter);
        scalar oneOverItP1 = 1./(avIter + 1);
        scalar mult = avIter*oneOverItP1;
        pMeanPtr_() == pMeanPtr_()*mult + pInst()*oneOverItP1;
        UMeanPtr_() == UMeanPtr_()*mult + UInst()*oneOverItP1;
        phiMeanPtr_() == phiMeanPtr_()*mult + phiInst()*oneOverItP1;
        RASModelVariables_().computeMeanFields();
        ++iAverageIter;
    }
}


void incompressibleVars::computeMeanUnsteadyFields()
{
    if (solverControl_.doAverageTime())
    {
        Info<< "Averaging fields" << endl;
        const scalar dt = mesh_.time().deltaTValue();
        const scalar elapsedTime
            = mesh_.time().value() - solverControl_.averageStartTime();
        scalar oneOverItP1 = dt/(elapsedTime + dt);
        scalar mult = elapsedTime/(elapsedTime + dt);
        pMeanPtr_() == pMeanPtr_()*mult + pInst()*oneOverItP1;
        UMeanPtr_() == UMeanPtr_()*mult + UInst()*oneOverItP1;
        phiMeanPtr_() == phiMeanPtr_()*mult + phiInst()*oneOverItP1;
        RASModelVariables_().computeMeanUnsteadyFields();
    }
}


void incompressibleVars::correctBoundaryConditions()
{
    correctNonTurbulentBoundaryConditions();
    correctTurbulentBoundaryConditions();
}


bool incompressibleVars::storeInitValues() const
{
    return solverControl_.storeInitValues();
}


bool incompressibleVars::computeMeanFields() const
{
    return solverControl_.average();
}


void incompressibleVars::transfer(variablesSet& vars)
{
    incompressibleVars& incoVars = refCast<incompressibleVars>(vars);
    // Copy source fields to the ones known by the object
    swapAndRename(pPtr_, incoVars.pPtr_);
    swapAndRename(UPtr_, incoVars.UPtr_);
    swapAndRename(phiPtr_, incoVars.phiPtr_);

    // Transfer turbulent fields. Copies fields since original fields are
    // not owned by RASModelVariables but from the turbulence model
    RASModelVariables_->transfer(incoVars.RASModelVariables()());
}


void incompressibleVars::validateTurbulence()
{
    // Update nut
    if (RASModelVariables_().hasNut())
    {
        turbulence()->validate();
    }
}


bool incompressibleVars::write() const
{
    // Write dummy fields, for continuation only
    if (useSolverNameForFields_)
    {
        if (RASModelVariables_().hasTMVar1())
        {
            createZeroFieldPtr<scalar>
            (
                mesh_,
                RASModelVariables_().TMVar1BaseName(),
                RASModelVariables_().TMVar1Inst().dimensions()
            )().write();
        }
        if (RASModelVariables_().hasTMVar2())
        {
            createZeroFieldPtr<scalar>
            (
                mesh_,
                RASModelVariables_().TMVar2BaseName(),
                RASModelVariables_().TMVar2Inst().dimensions()
            )().write();
        }
        if (RASModelVariables_().hasNut())
        {
            createZeroFieldPtr<scalar>
            (
                mesh_,
                RASModelVariables_().nutBaseName(),
                RASModelVariables_().nutRefInst().dimensions()
            )().write();
        }

        return true;
    }

    return false;
}


fvMesh& incompressibleVars::mesh() const
{
    return mesh_;
}


const solverControl& incompressibleVars::solverControlReference() const
{
    return solverControl_;
}


const wordList& incompressibleVars::allocatedFieldNames() const
{
    return allocatedFieldNames_;
}


const Vector<label>& incompressibleVars::solDirs() const
{
    return solDirs_;
}


void incompressibleVars::adjustAverageStartTime(const scalar& offset)
{
    solverControl_.averageStartTime() += offset;
}


void incompressibleVars::setWriteOption(IOobject::writeOption w)
{
    if (!writeFields_)
    {
        w = IOobject::NO_WRITE;
    }
    pInst().writeOpt() = w;
    phiInst().writeOpt() = w;
    UInst().writeOpt() = w;
    if (solverControl_.doAverageTime())
    {
        pMeanPtr_->writeOpt() = w;
        phiMeanPtr_->writeOpt() = w;
        UMeanPtr_->writeOpt() = w;
    }
    incompressible::RASModelVariables& rasVars = RASModelVariables_();
    rasVars.setWriteOption(w);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
