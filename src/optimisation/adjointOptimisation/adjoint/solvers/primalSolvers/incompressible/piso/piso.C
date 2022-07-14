/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
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

#include "piso.H"
#include "findRefCell.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "Time.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"
#include "profiling.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(piso, 0);
    addToRunTimeSelectionTable(incompressiblePrimalSolver, piso, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::piso::printFields()
{
    if (debug > 1)
    {
        const word timeName = mesh_.time().timeName();
        const word prefix = "after solving for t = " + timeName + " ";
        volScalarField& p = incoVars_.pInst();
        volVectorField& U = incoVars_.UInst();
        surfaceScalarField& phi = incoVars_.phiInst();
        Info<< "Debug: U.nOldTimes()   = " << U.nOldTimes() << endl;
        Info<< "Debug: phi.nOldTimes() = " << phi.nOldTimes() << endl;
        Info<< prefix << "p:\n" << p << endl;
        Info<< prefix << "U:\n" << U << endl;
        Info<< prefix << "U.oldTime():\n" << U.oldTime() << endl;
        Info<< prefix << "phi:\n" << phi << endl;
        Info<< prefix << "phi.oldTime():\n" << phi.oldTime() << endl;
        incompressible::RASModelVariables& rasVars =
            incoVars_.RASModelVariables()();
        if (rasVars.hasTMVar1())
        {
            Info<< "Debug: TMVar1.nOldTimes() = "
                << rasVars.TMVar1Inst().nOldTimes() << endl;
            Info<< prefix << " TMVar1:\n" << rasVars.TMVar1Inst() << endl;
            Info<< prefix << " TMVar1.oldTime():\n"
                << rasVars.TMVar1Inst().oldTime() << endl;
        }
        if (rasVars.hasTMVar2())
        {
            Info<< "Debug: TMVar2.nOldTimes() = "
                << rasVars.TMVar2Inst().nOldTimes() << endl;
            Info<< prefix << " TMVar2:\n" << rasVars.TMVar2Inst() << endl;
            Info<< prefix << " TMVar2.oldTime():\n"
                << rasVars.TMVar2Inst().oldTime() << endl;
        }
        if (rasVars.hasNut())
        {
            Info << prefix << " nut:\n" << rasVars.nutRefInst() << endl;
        }
    }
}


Foam::incompressibleVars& Foam::piso::allocateVars(scalar readTime)
{
    vars_.reset(new incompressibleVars(mesh_, solverControl_(), readTime));
    return getIncoVars();
}


void Foam::piso::CourantNo() const
{
    const surfaceScalarField& phi = incoVars_.phiInst();
    const Time& runTime = mesh_.time();
    const fvMesh& mesh = mesh_;
    #include "CourantNo.H"
}


void Foam::piso::continuityErrors()
{
    // Original function at
    // src/finiteVolume/cfdTools/incompressible/continuityErrs.H
    const Time& runTime = mesh_.time();
    const fvMesh& mesh = mesh_;
    const surfaceScalarField& phi = incoVars_.phiInst();
    scalar& cumulativeContErr = cumulativeContErr_;
    #include "continuityErrs.H"
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::piso::piso
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    Switch useCustomReadTime
)
:
    incompressiblePrimalSolver
    (
        mesh,
        managerType,
        dict,
        useCustomReadTime
    ),
    solverControl_(PISOControl::New(mesh, managerType, *this)),
    shouldResetMeanFields_(useCustomReadTime),
    incoVars_
    (
        allocateVars
        (
            useCustomReadTime ? dict.getOrDefault<scalar>("readTime", 0) : -1
        )
    ),
    MRF_(mesh, word(useSolverNameForFields() ? solverName() : word::null)),
    cumulativeContErr_(Zero),
    objectives_(0),
    timeManip_
    (
        const_cast<unsteadyTimeManipulation&>
        (
            unsteadyTimeManipulation::New(mesh)
        )
    )
{
    addExtraSchemes();
    setRefCell
    (
        incoVars_.pInst(),
        solverControl_().dict(),
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );

    // Decide primal storage mechanism
    if (dict.isDict("storage"))
    {
        primalStorage_.reset
        (
            primalStorage::New(*this, dict.subDict("storage"))
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::piso::readDict(const dictionary& dict)
{
    if (incompressiblePrimalSolver::readDict(dict))
    {
        return true;
    }

    return false;
}


void Foam::piso::solveIter()
{
    addProfiling(piso, "piso::solveIter()");
    const Time& time = mesh_.time();
    // The primal solution at the first time-step does not need to be stored.
    // However, when checkpointing is used, the checkpoint at the first time
    // step must be set; storeInitialVariables() does nothing in other cases.
    // Note: would make sense to be at the end of the iteration but
    // checkpoint intricaties suggest it should remain here
    if (primalStorage_.valid())
    {
        if (time.timeIndex() == timeManip_.startTimeIndex() + 1)
        {
            primalStorage_().storeInitialVariables();
        }
        else
        {
            primalStorage_().storeVariables();
        }
    }
    Info<< "Time = " << mesh_.time().timeName() << "\n" << endl;

    // Grab references
    volScalarField& p = incoVars_.pInst();
    volVectorField& U = incoVars_.UInst();
    surfaceScalarField& phi = incoVars_.phiInst();
    autoPtr<incompressible::turbulenceModel>& turbulence =
        incoVars_.turbulence();
    label&  pRefCell = solverControl_().pRefCell();
    scalar& pRefValue = solverControl_().pRefValue();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    CourantNo();

    // Momentum predictor
    //~~~~~~~~~~~~~~~~~~~

    MRF_.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF_.DDt(U)
      + turbulence->divDevReff(U)
      ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (solverControl_().momentumPredictor())
    {
        Foam::solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }

    // Pressure Eq
    //~~~~~~~~~~~~
    while (solverControl_().correct())
    {
        volScalarField rAU(1.0/UEqn.A());
        volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::flux(HbyA)
          + MRF_.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi))
        );

        MRF_.makeRelative(phiHbyA);

        adjustPhi(phiHbyA, U, p);

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p, U, phiHbyA, rAU, MRF_);

        // Non-orthogonal pressure corrector loop
        while (solverControl_().correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
            );

            pEqn.setReference(pRefCell, pRefValue);

            fvOptions.constrain(pEqn);

            pEqn.solve
            (
                mesh_.solver(p.select(solverControl_().finalInnerIter()))
            );

            if (solverControl_().finalNonOrthogonalIter())
            {
                phi = phiHbyA - pEqn.flux();
            }
        }

        continuityErrors();

        U = HbyA - rAU*fvc::grad(p);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
    // An extra correction of the BC of p would be necessary to have
    // consistency for p upon reconstruction during the adjoint run.
    // However, this would slightly affect the primal solution, which is
    // something we want to avoid. A small inconsistency between the field
    // computed by the primal solver and the reconstructed one may occur, when
    // BCs depending on the flux value (e.g.  outletInlet) are used. To avoid
    // this inconsistency activate 'storeAllBoundaries' in the storade
    // subdictionary, although this means that all field boundaries will be
    // stored and the correctBCs will not be called upon reconstruction.
    // p.correctBoundaryConditions();

    incoVars_.laminarTransport().correct();
    turbulence->correct();

    // Print objective values to screen and compute mean value
    Info<< endl;
    for (objective* obj : objectives_)
    {
        Info<< obj->objectiveName() << " : " << obj->J() << endl;
        obj->accumulateJMean();
        obj->writeInstantaneousValue();
    }

    // Average fields if necessary
    incoVars_.computeMeanUnsteadyFields();

    // Calls time.write()
    solverControl_().writeFields();

    // Print execution time
    time.printExecutionTime(Info);

    printFields();
}


void Foam::piso::solve()
{
    // Iterate
    if (active_)
    {
        preLoop();
        while (loop())
        {
            this->solveIter();
        }
        postLoop();
    }
}


bool Foam::piso::loop()
{
    Time& time = const_cast<Time&>(mesh_.time());
    return time.loop();
}


void Foam::piso::restoreInitValues()
{
    incoVars_.restoreInitValues();
}


void Foam::piso::preLoop()
{
    // Move to the correct time, if changed by a previous primal solver
    timeManip_.moveToPrimalStartTime();

    // Get the objectives for this solver
    if (objectives_.empty())
    {
        objectives_ = getObjectiveFunctions();
    }

    // Reset initial and mean fields before solving
    restoreInitValues();
    if (shouldResetMeanFields_)
    {
        incoVars_.resetMeanFields();
    }

    // Validate turbulence model fields
    incoVars_.turbulence()->validate();

    // nullify JMean for all objectives
    // Note: we should cover the same of continuing avegaring the objective
    // after a primal restart
    for (objective* obj : objectives_)
    {
        obj->resetJMean();
    }

    // Reset primal storage
    if (primalStorage_.valid())
    {
        primalStorage_().scratch();
    }

    // Print out
    Info<< nl << "- - - - - - - - - - - - - - - - - - - - -" << endl;
    Info<< "Solving primal equations for solver " << solverName() << endl;
    Info<< "- - - - - - - - - - - - - - - - - - - - -" << nl << endl;

    printFields();
}


void Foam::piso::postLoop()
{
    addProfiling(piso, "piso::postLoop()");
    // Store last time instance. Redundant if a single adjoint solver is used,
    // but necessary for more than one adjoint solvers
    if (primalStorage_.valid())
    {
        // Advance time to avoid implications with the checkpointing
        // implementation and store variables
        timeManip_.storeTime();
        const_cast<Time&>(mesh_.time())++;
        primalStorage_().storeVariables();

        // Restore time
        timeManip_.restoreTime();

        primalStorage_().storageMetrics();
    }
    // Write separators in objective files
    for (objective* obj : objectives_)
    {
        obj->writeInstantaneousSeparator();
    }
    // Safety
    objectives_.clear();

    // Avoid writting fields if already in an outputTime iter
    // since results will be written by the solver class either way
    Time& time = const_cast<Time&>(mesh_.time());
    if (!time.writeTime())
    {
        // Write all fields
        time.writeNow();
        // Write dummy turbulence fields if multiple primal solvers exist,
        // to support seamless restarts
        this->writeNow();
    }

    // Allow mean fields to be zeroed at the next call of preLoop
    shouldResetMeanFields_ = true;

    // Store endTimeIndex, to be retrieved later by the adjoint
    timeManip_.storeEndTimeIndex();
}


void Foam::piso::postLineSearch()
{
    // Update averaging window for the next optimisation cycle
    solverControl_().updateAverageStartTime();
}


bool Foam::piso::writeData(Ostream& os) const
{
    os.writeEntry("averageIter", solverControl_().averageIter());

    return true;
}


// ************************************************************************* //
