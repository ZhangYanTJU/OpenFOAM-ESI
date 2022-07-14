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

#include "pimple.H"
#include "findRefCell.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "CorrectPhi.H"
#include "Time.H"
#include "fvOptions.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pimple, 0);
    addToRunTimeSelectionTable
    (
        incompressiblePrimalSolver,
        pimple,
        dictionary
    );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::incompressibleVars& Foam::pimple::allocateVars(scalar readTime)
{
    vars_.reset(new incompressibleVars(mesh_, solverControl_(), readTime));
    return getIncoVars();
}


void Foam::pimple::readTimeControls()
{
    // Original function at
    // src/finiteVolume/cfdTools/general/include/readTimeControls.H
    const Time& runTime = mesh_.time();
    bool& adjustTimeStep = adjustTimeStep_;
    scalar& maxCo = maxCo_;
    scalar& maxDeltaT = maxDeltaT_;
    #include "readTimeControls.H"
}


void Foam::pimple::CourantNo()
{
    // Original function at
    // src/finiteVolume/cfdTools/incompressible/CourantNo.H
    const surfaceScalarField& phi = incoVars_.phiInst();
    const Time& runTime = mesh_.time();
    const fvMesh& mesh = mesh_;
    #include "CourantNo.H"
    CoNum_ = CoNum;
}


void Foam::pimple::continuityErrors()
{
    // Original function at
    // src/finiteVolume/cfdTools/incompressible/continuityErrs.H
    const Time& runTime = mesh_.time();
    const fvMesh& mesh = mesh_;
    const surfaceScalarField& phi = incoVars_.phiInst();
    scalar& cumulativeContErr = cumulativeContErr_;
    #include "continuityErrs.H"
}


void Foam::pimple::setDeltaT()
{
    // Original function at
    // src/finiteVolume/cfdTools/general/include/setDeltaT.H
    Time& runTime = const_cast<Time&>(mesh_.time());
    bool& adjustTimeStep = adjustTimeStep_;
    scalar& maxCo = maxCo_;
    scalar& maxDeltaT = maxDeltaT_;
    scalar& CoNum = CoNum_;
    #include "setDeltaT.H"
}


void Foam::pimple::setInitialDeltaT()
{
    // Original function at
    // src/finiteVolume/cfdTools/general/include/setInitialDeltaT.H
    Time& runTime = const_cast<Time&>(mesh_.time());
    const bool& adjustTimeStep = adjustTimeStep_;
    const scalar& CoNum = CoNum_;
    const scalar& maxCo = maxCo_;
    const scalar& maxDeltaT = maxDeltaT_;
    #include "setInitialDeltaT.H"
}


void Foam::pimple::setTime()
{
    Time& runTime = const_cast<Time&>(mesh_.time());
    ++runTime;
}


void Foam::pimple::storeDeltaT()
{
    Time& time = const_cast<Time&>(mesh_.time());
    if (adjustTimeStep_)
    {
        DebugInfo
            << "time.timeIndex() - primalStartTimeIndex = "
            << time.timeIndex() - timeManip_.startTimeIndex() << nl
            << "time.timeIndex(), primalStartTimeIndex  = "
            << time.timeIndex() << " " << timeManip_.startTimeIndex()
            << endl;
        if
        (
            deltaTList_.empty()
         && (time.timeIndex() > timeManip_.startTimeIndex())
        )
        {
            deltaTList_.
                setSize
                (
                    time.timeIndex() - timeManip_.startTimeIndex(),
                    deltaT0_
                );
        }
        deltaTList_.append(time.deltaTValue());
    }
    DebugInfo
        << "deltaTList.size() = " << deltaTList_.size() << endl;
}


void Foam::pimple::scratchDeltaTHistory()
{
    const Time& time = mesh_.time();
    deltaT0_ = time.deltaTValue();
    primalStartTimeIndex_ = time.timeIndex();
    deltaTList_.clear();
}


const Foam::DynamicList<Foam::scalar>& Foam::pimple::getDeltaTList() const
{
    return deltaTList_;
}


/*void Foam::pimple::readDyMControls()
{
    // Original function at src/dynamicFvMesh/include/readDyMControls.H
    // Original function at
    // src/finiteVolume/cfdTools/general/include/readTimeControls.H
    const Time& runTime = mesh_.time();
    const PIMPLEControl& pimple = solverControl_();
    bool& moveMeshOuterCorrectors = moveMeshOuterCorrectors_;
    bool& checkMeshCourantNo = checkMeshCourantNo_;
    bool& correctPhi = correctPhi_;
    bool& adjustTimeStep = adjustTimeStep_;
    scalar& maxCo = maxCo_;
    scalar& maxDeltaT = maxDeltaT_;
    #include "readDyMControls.H"
}


void Foam::pimple::meshCourantNo()
{
    // Original function at src/dynamicFvMesh/include/meshCourantNo.H
    const Time& runTime = mesh_.time();
    const fvMesh& mesh = mesh_;
    #include "meshCourantNo.H"
}


void Foam::pimple::correctPhi()
{
    // Original function at the location of the solver
    // CorrectPhi at src/finiteVolume/cfdTools/general/CorrectPhi/CorrectPhi.H
    volScalarField& p = incoVars_.pInst();
    volVectorField& U = incoVars_.UInst();
    surfaceScalarField& phi = incoVars_.phiInst();
    CorrectPhi
    (
        U,
        phi,
        p,
        dimensionedScalar("rAUf", dimTime, 1),
        geometricZeroField(),
        solverControl_()
    );
    continuityErrors();
}


void Foam::pimple::createUfIfPresent()
{
    // Original function at
    // src/finiteVolume/cfdTools/incompressible/createUfIfPresent.H
    if (mesh_.dynamic())
    {
        volVectorField& U = incoVars_.UInst();

        Info<< "Constructing face velocity Uf\n" << endl;

        Uf_.reset
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "Uf",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U)
            )
        );
    }
}*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pimple::pimple
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
    solverControl_(PIMPLEControl::New(mesh, managerType, *this)),
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
    adjustTimeStep_
    (
        mesh.time().controlDict().getOrDefault("adjustTimeStep", false)
    ),
    maxCo_
    (
        mesh.time().controlDict().getOrDefault<scalar>("maxCo", 1.0)
    ),
    maxDeltaT_
    (
        mesh.time().controlDict().getOrDefault<scalar>("maxDeltaT", GREAT)
    ),
    primalStartTimeIndex_(mesh_.time().timeIndex()),
    deltaTList_(0),
    timeManip_
    (
        const_cast<unsteadyTimeManipulation&>
        (
            unsteadyTimeManipulation::New(mesh)
        )
    )
    /*correctPhi_
    (
        solverControl_().dict().getOrDefault("correctPhi", mesh_.dynamic())
    ),
    checkMeshCourantNo_
    (
        solverControl_().dict().getOrDefault("checkMeshCourantNo", false)
    ),
    moveMeshOuterCorrectors_
    (
        solverControl_().dict().getOrDefault("moveMeshOuterCorrectors", false)
    )*/
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

bool Foam::pimple::readDict(const dictionary& dict)
{
    if (incompressiblePrimalSolver::readDict(dict))
    {
        return true;
    }

    return false;
}


void Foam::pimple::solveIter()
{
    const Time& time = mesh_.time();
    //readDyMControls();
    readTimeControls();
    setDeltaT();
    storeDeltaT();
    setTime();

    // The primal solution at the first time-step does not need to be stored.
    // When checkpointing is used, the checkpoint at the first time step must
    // be set, else storeInitialVariables() does nothing
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

    while (solverControl_().loop())
    {
        // No moving geometries supported for now
        /*
        if (solverControl_().firstIter() || moveMeshOuterCorrectors_)
        {
            //mesh_.update();

            if (mesh_.changing())
            {
                MRF_.update();

                if (correctPhi_)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                    phi = mesh_.Sf() & Uf_();

                    correctPhi();

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi, U);
                }
                if (checkMeshCourantNo_)
                {
                    meshCourantNo();
                }
            }
        }
        */

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
              + MRF_.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf_))
            );

            MRF_.makeRelative(phiHbyA);

            if (p.needReference())
            {
                fvc::makeRelative(phiHbyA, U);
                adjustPhi(phiHbyA, U, p);
                fvc::makeAbsolute(phiHbyA, U);
            }

            tmp<volScalarField> rAtU(rAU);

            if (solverControl_().consistent())
            {
                rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
                phiHbyA +=
                    fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh_.magSf();
                HbyA -= (rAU - rAtU())*fvc::grad(p);
            }

            if (solverControl_().nCorrPISO() <= 1)
            {
                tUEqn.clear();
            }

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAtU(), MRF_);

            // Non-orthogonal pressure corrector loop
            while (solverControl_().correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
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

            // Explicitly relax pressure for momentum corrector
            p.relax();

            U = HbyA - rAtU*fvc::grad(p);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            // Correct Uf if the mesh is moving
            fvc::correctUf(Uf_, U, phi);

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);
        }

        if (solverControl_().turbCorr())
        {
            incoVars_.laminarTransport().correct();
            turbulence->correct();
        }
    }

    // Print objective values to screen and compute mean value
    Info<< nl;
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
}


void Foam::pimple::solve()
{
    // Iterate
    if (active_)
    {
        preLoop();
        while (loop())
        {
            solveIter();
        }
        postLoop();
    }
}


bool Foam::pimple::run()
{
  //return solverControl_().run();
    return mesh_.time().run();
}


bool Foam::pimple::loop()
{
    // mesh_.time().run() is actually needed here. Not a mistake.
    // Chosen for compatibility reasons with the simple and piso algorithms
    // and the 'solveWithArgs' function at solverTemplates.C
    return mesh_.time().run();
}


void Foam::pimple::restoreInitValues()
{
    incoVars_.restoreInitValues();
}


void Foam::pimple::preLoop()
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
    for (objective* obj : objectives_)
    {
        obj->resetJMean();
    }

    // Reset primal storage
    if (primalStorage_.valid())
    {
        primalStorage_().scratch();
    }
    scratchDeltaTHistory();
    CourantNo();
    //createUfIfPresent();
    setInitialDeltaT();

    // Print out
    Info<< nl << "- - - - - - - - - - - - - - - - - - - - -" << endl;
    Info<< "Solving primal equations for solver " << solverName() << endl;
    Info<< "- - - - - - - - - - - - - - - - - - - - -" << nl << endl;
}


void Foam::pimple::postLoop()
{
    // Store last time instance. Redundant if a single adjoint solver is used,
    // but necessary for more than one adjoint solvers
    if (primalStorage_.valid())
    {
        // Advance time to avoid implications with the checkpointing
        // implementation and store variables
        timeManip_.storeTime();
        const_cast<Time&>(mesh_.time())++;
        primalStorage_().storeVariables();

        // Retrive time
        timeManip_.restoreTime();

        // Print storage metrics
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
        // Write dummy turbulence fields if multiple primal solvers exist
        this->writeNow();
    }

    // Allow mean fields to be zeroed at the next call of preLoop
    shouldResetMeanFields_ = true;

    // Store endTimeIndex, to be retrieved later by the adjoint
    timeManip_.storeEndTimeIndex();
}


void Foam::pimple::postLineSearch()
{
    // Update averaging window for the next optimisation cycle
    solverControl_().updateAverageStartTime();
}


bool Foam::pimple::writeData(Ostream& os) const
{
    os.writeEntry("averageIter", solverControl_().averageIter());

    return true;
}

// ************************************************************************* //
