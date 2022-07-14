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

#include "pimple.H"
#include "adjointPimple.H"
#include "findRefCell.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(adjointPimple, 0);
    addToRunTimeSelectionTable
    (
        incompressibleAdjointSolver,
        adjointPimple,
        dictionary
    );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::incompressibleAdjointVars& Foam::adjointPimple::allocateVars()
{
    vars_.reset
    (
        new incompressibleAdjointVars
        (
            mesh_,
            solverControl_(),
            objectiveManagerPtr_(),
            primalVars_
        )
    );
    return getAdjointVars();
}


void Foam::adjointPimple::setDeltaT()
{
    // Original function at
    // src/finiteVolume/cfdTools/general/include/setDeltaT.H
    if (!deltaTList_.empty())
    {
        Time& runTime = const_cast<Time&>(mesh_.time());
        DebugInfo
            << "runTime.timeIndex(), primalStartTimeIndex = "
            << runTime.timeIndex() << " " << timeManip_.startTimeIndex()
            << endl;
        runTime.setDeltaT
        (
            deltaTList_[runTime.timeIndex() - timeManip_.startTimeIndex() -1]
        );

        Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
    }

}


void Foam::adjointPimple::setTime()
{
    Time& time = const_cast<Time&>(mesh_.time());
    --time;
}


void Foam::adjointPimple::CourantNo() const
{
    const surfaceScalarField& phi = primalVars_.phiInst();
    const Time& runTime = mesh_.time();
    const fvMesh& mesh = mesh_;
    #include "CourantNo.H"
}


void Foam::adjointPimple::continuityErrors()
{
    // Original function at
    // src/finiteVolume/cfdTools/incompressible/continuityErrs.H
    const Time& runTime = mesh_.time();
    const fvMesh& mesh = mesh_;
    const surfaceScalarField& phi = adjointVars_.phiaInst();
    scalar& cumulativeContErr = cumulativeContErr_;
    #include "continuityErrs.H"
}


void Foam::adjointPimple::accumulateSensitivities()
{
    if (computeSensitivities_)
    {
        adjointSensitivity_->accumulateIntegrand(mesh_.time().deltaTValue());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointPimple::adjointPimple
(
    fvMesh& mesh,
    const word& managerType,
    const dictionary& dict,
    const word& primalSolverName
)
:
    incompressibleAdjointSolver
    (
        mesh,
        managerType,
        dict,
        primalSolverName
    ),
    solverControl_(PIMPLEControl::New(mesh, managerType, *this)),
    adjointVars_(allocateVars()),
    cumulativeContErr_(Zero),
    adjointSensitivity_(nullptr),
    primalStorage_(getPrimalSolver().getPrimalStorage()),
    deltaTList_(refCast<pimple>(getPrimalSolver()).getDeltaTList()),
    timeManip_
    (
        const_cast<unsteadyTimeManipulation&>
        (
            unsteadyTimeManipulation::New(mesh)
        )
    )
{
    ATCModel_.reset
    (
        ATCModel::New
        (
            mesh,
            primalVars_,
            adjointVars_,
            dict.subDict("ATCModel")
        ).ptr()
    );

    addExtraSchemes();
    setRefCell
    (
        adjointVars_.paInst(),
        solverControl_().dict(),
        solverControl_().pRefCell(),
        solverControl_().pRefValue()
    );

    if (computeSensitivities_)
    {
        const IOdictionary& optDict =
            mesh.lookupObject<IOdictionary>("optimisationDict");

        adjointSensitivity_.reset
        (
            incompressible::adjointSensitivity::New
            (
                mesh,
                optDict.subDict("optimisation").subDict("sensitivities"),
                *this
            ).ptr()
        );
    }

    // Check integration times of the objectives
    objectiveManagerPtr_().checkIntegrationTimes();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adjointPimple::readDict(const dictionary& dict)
{
    if (incompressibleAdjointSolver::readDict(dict))
    {
        if (adjointSensitivity_.valid())
        {
            const IOdictionary& optDict =
                mesh_.lookupObject<IOdictionary>("optimisationDict");

            adjointSensitivity_().readDict
            (
                optDict.subDict("optimisation").subDict("sensitivities")
            );
        }

        return true;
    }

    return false;
}


void Foam::adjointPimple::solveIter()
{
    CourantNo();
    setDeltaT();
    setTime();

    // Retrieve primal fields of the time step. Retrieving of the first
    // time-step of the adjoint loop is also necessary if more than one adjoint
    // solvers are used. Includes potential recomputation of the primal fields,
    // depending on the primal storage strategy
    primalStorage_->retrieveVariables();
    const Time& time = mesh_.time();
    Info<< "Time = " << time.timeName() << "\n" << endl;

    // Update objective related fields
    objectiveManagerPtr_->updateOrNullify();

    // Update turbulence model and ATC related terms
    incompressibleAdjointSolver::updatePrimalBasedQuantities();

    // Grab primal references
    const surfaceScalarField& phi = primalVars_.phi();
    // Grab adjoint references
    volScalarField& pa = adjointVars_.paInst();
    volVectorField& Ua = adjointVars_.UaInst();
    surfaceScalarField& phia = adjointVars_.phiaInst();
    autoPtr<incompressibleAdjoint::adjointRASModel>& adjointTurbulence =
        adjointVars_.adjointTurbulence();
    const label&  paRefCell  = solverControl_().pRefCell();
    const scalar& paRefValue = solverControl_().pRefValue();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    while (solverControl_().loop())
    {
        // Momentum predictor
        //~~~~~~~~~~~~~~~~~~~

        tmp<fvVectorMatrix> tUaEqn
        (
            fvm::ddt(Ua) + fvm::div(-phi, Ua)
          + adjointTurbulence->divDevReff(Ua)
          + adjointTurbulence->adjointMeanFlowSource()
          ==
            fvOptions(Ua)
        );
        fvVectorMatrix& UaEqn = tUaEqn.ref();

        // Add sources from boundary conditions
        UaEqn.boundaryManipulate(Ua.boundaryFieldRef());

        // Add sources from volume-based objectives
        objectiveManagerPtr_().addUaEqnSource(UaEqn);

        // Add ATC term
        ATCModel_->addATC(UaEqn);

        UaEqn.relax();

        fvOptions.constrain(UaEqn);

        if (solverControl_().momentumPredictor())
        {
            Foam::solve(UaEqn == -fvc::grad(pa));

            fvOptions.correct(Ua);
        }

        // Pressure Eq
        //~~~~~~~~~~~~
        while (solverControl_().correct())
        {
            volScalarField rAUa(1.0/UaEqn.A());
            volVectorField HabyA(constrainHbyA(rAUa*UaEqn.H(), Ua, pa));
            surfaceScalarField phiaHbyA
            (
                "phiaHbyA",
                fvc::flux(HabyA)
              //+ MRF_.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi))
            );

            if (pa.needReference())
            {
                fvc::makeRelative(phiaHbyA, Ua);
                adjustPhi(phiaHbyA, Ua, pa);
                fvc::makeAbsolute(phiaHbyA, Ua);
            }

            tmp<volScalarField> rAtUa(rAUa);

            if (solverControl_().consistent())
            {
                rAtUa = 1.0/max(1.0/rAUa - UaEqn.H1(), 0.1/rAUa);
                phiaHbyA +=
                    fvc::interpolate(rAtUa() - rAUa)*fvc::snGrad(pa)*mesh_.magSf();
                HabyA -= (rAUa - rAtUa())*fvc::grad(pa);
            }

            if (solverControl_().nCorrPISO() <= 1)
            {
                tUaEqn.clear();
            }

            // Non-orthogonal pressure corrector loop
            while (solverControl_().correctNonOrthogonal())
            {
                fvScalarMatrix paEqn
                (
                    fvm::laplacian(rAtUa(), pa) == fvc::div(phiaHbyA)
                );

                paEqn.boundaryManipulate(pa.boundaryFieldRef());

                fvOptions.constrain(paEqn);
                paEqn.setReference(paRefCell, paRefValue);

                paEqn.solve
                (
                    mesh_.solver(pa.select(solverControl_().finalInnerIter()))
                );

                if (solverControl_().finalNonOrthogonalIter())
                {
                    phia = phiaHbyA - paEqn.flux();
                }
            }

            continuityErrors();

            // Explicitly relax pressure for momentum corrector
            pa.relax();

            Ua = HabyA - rAtUa()*fvc::grad(pa);
            Ua.correctBoundaryConditions();
            fvOptions.correct(Ua);
            pa.correctBoundaryConditions();
        }

        if (solverControl_().turbCorr())
        {
            adjointTurbulence->correct();
        }

        if (solverControl_().printMaxMags())
        {
            dimensionedScalar maxUa = max(mag(Ua)());
            dimensionedScalar maxpa = max(mag(pa)());
            Info<< "Max mag of adjoint velocity = " << maxUa.value() << endl;
            Info<< "Max mag of adjoint pressure = " << maxpa.value() << endl;
        }
    }

    // Average fields if necessary
    adjointVars_.computeMeanUnsteadyFields(timeManip_.primalEndTime());

    // Calls time.write()
    solverControl_().writeFields();

    // Print execution time
    time.printExecutionTime(Info);
}


void Foam::adjointPimple::solve()
{
    // Iterate
    if (active_)
    {
        preLoop();
        while (loop())
        {
            solveIter();
            accumulateSensitivities();
        }
        postLoop();
    }
}


bool Foam::adjointPimple::run()
{
    Time& time = const_cast<Time&>(mesh_.time());
    return time.reverseRun();
}


bool Foam::adjointPimple::loop()
{
    Time& time = const_cast<Time&>(mesh_.time());
    // mesh_.time().reverseRun() is actually needed here. Not a mistake.
    // Chosen for compatibility reasons with the piso algorithm
    // and the 'solveWithArgs' function at solverTemplates.C
    return time.reverseRun();
}


void Foam::adjointPimple::preLoop()
{
    // Disable re-writting of primal fields
    primalVars_.setWriteOption(IOobject::NO_WRITE);

    // Reset time to the beginning of the adjoint loop.
    // Resets endTime too
    timeManip_.moveToAdjointStartTime();

    //  TODO: To be able to stop the adjoint solution midway and continue
    //  thereafter, the oldTimes of the adjoint fields are required to be
    //  correctly read.  Make sure that these are not reset due to change in
    //  time, upon calling oldTime() function.

    //  Reset all adjoint fields to zero
    adjointVars_.nullify();
    adjointVars_.resetMeanFields();

    // Copy pointers to storage instances and possibly read compressed fields
    // in cases of restarted runs
    primalStorage_->preAdjointLoop();

    // The adjoint equations should be also solved at the last primal time-step.
    // Advance time once more to get to the correct time during the solution
    // of the first adjoint time-setp
    const_cast<Time&>(mesh_.time())++;

    // Print out
    Info<< nl << "- - - - - - - - - - - - - - - - - - - - -" << endl;
    Info<< "Solving adjoint equations for solver " << solverName() << endl;
    Info<< "- - - - - - - - - - - - - - - - - - - - -" << nl << endl;
}


void Foam::adjointPimple::postLoop()
{
    // rite adjoint fields at primalEndTime. Although the adjoint equations
    // are integrated backwards in time starting from primalEndTime and ending
    // at primalStartTime, the appropriate time to write the adjoint fields is
    // the primalEndTime. This is not only to have both the primal and adjoint
    // fields of each optimization cycle written at the same time, but more
    // importantly to have consistency between the written grid and the adjoint
    // fields. To do so, reset only the time value to that corresponding to the
    // end of the primal. Keep the timeIndex unchanged, so that
    // GeometricField::storeOldTime() function is not called by the time
    // GeometricField::oldTime() is called from the
    // incompressibleAdjointVars::write() function. This way, the adjoint
    // fields will be written correctly in disk without oldTimes to been
    // shifted
    Time& time = const_cast<Time&>(mesh_.time());
    timeManip_.moveToAdjointStartTime(false, false);
    if (!time.writeTime())
    {
        adjointVars_.write();
    }
    // Completely reset time to that corresponding to the end of the primal
    // problem.
    timeManip_.moveToAdjointStartTime();
    // Write data except for the adjoint fields that have already been written
    if (!time.writeTime())
    {
        adjointVars_.setWriteOption(IOobject::NO_WRITE);
        time.writeNow();
        this->writeNow();
        adjointVars_.setWriteOption(IOobject::AUTO_WRITE);
    }

    // Print useful info and rewind storage pointers
    primalStorage_->postAdjointLoop();

    // Compute sensitivities
    adjointSolver::postLoop();

    //- Re-enable writting of primal fields
    primalVars_.setWriteOption(IOobject::AUTO_WRITE);
}


void Foam::adjointPimple::postLineSearch()
{
    // Set integration times for the next optimisation cycle
    objectiveManagerPtr_().incrementIntegrationTimes(timeManip_.span());
}


void Foam::adjointPimple::computeObjectiveSensitivities()
{
    if (computeSensitivities_)
    {
        const scalarField& sens = adjointSensitivity_->calculateSensitivities();
        if (!sensitivities_)
        {
            sensitivities_.reset(new scalarField(sens.size(), scalar(0)));
        }
        sensitivities_.ref() = sens;
    }
    else
    {
        sensitivities_.reset(new scalarField(0));
    }
}


const Foam::scalarField& Foam::adjointPimple::getObjectiveSensitivities()
{
    if (!sensitivities_.valid())
    {
        computeObjectiveSensitivities();
    }

    return sensitivities_();
}


void Foam::adjointPimple::clearSensitivities()
{
    if (computeSensitivities_)
    {
        adjointSensitivity_->clearSensitivities();
        adjointSolver::clearSensitivities();
    }
}


Foam::sensitivity& Foam::adjointPimple::getSensitivityBase()
{
    if (adjointSensitivity_.valid())
    {
        return adjointSensitivity_();
    }
    else
    {
        FatalErrorInFunction
            << "Sensitivity object not allocated \n"
            << "Turn computeSensitivities on in "
            << solverName_
            << nl << nl
            << exit(FatalError);

        return adjointSensitivity_();
    }
}


void Foam::adjointPimple::updatePrimalBasedQuantities()
{
    // Update turbulence model and ATC related terms
    incompressibleAdjointSolver::updatePrimalBasedQuantities();

    // Write mean values for all objectives
    PtrList<objective>& objectives
        = objectiveManagerPtr_->getObjectiveFunctions();
    for (objective& obj : objectives)
    {
        if (obj.shouldWrite())
        {
            obj.writeMeanValue();
        }
    }
}


bool Foam::adjointPimple::writeData(Ostream& os) const
{
    os.writeEntry("averageIter", solverControl_().averageIter());
    return adjointSolver::writeData(os);
}


// ************************************************************************* //
