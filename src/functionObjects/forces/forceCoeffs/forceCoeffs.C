/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "forceCoeffs.H"
#include "dictionary.H"
#include "Time.H"
#include "Pstream.H"
#include "IOmanip.H"
#include "fvMesh.H"
#include "dimensionedTypes.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forceCoeffs, 0);
    addToRunTimeSelectionTable(functionObject, forceCoeffs, dictionary);
}
}

Foam::HashTable<Foam::functionObjects::forceCoeffs::coeffDesc>&
Foam::functionObjects::forceCoeffs::coeffsSelection()
{
    static HashTable<coeffDesc> coeffs;

    if (coeffs.empty())
    {
        coeffs.insert("Cd", coeffDesc("Drag force coefficient", "Cd", 0, 0));
        coeffs.insert("Cs", coeffDesc("Side force coefficient", "Cs", 1, 2));
        coeffs.insert("Cl", coeffDesc("Lift force coefficient", "Cl", 2, 1));

        coeffs.insert("CmRoll", coeffDesc("Roll moment coefficient", "CmRoll", 0, -1));
        coeffs.insert("CmYaw", coeffDesc("Yaw moment coefficient", "CmYaw", 1, -1));
        coeffs.insert("CmPitch", coeffDesc("Pitch moment coefficient", "CmPitch", 2, -1));
    }

    return coeffs;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forceCoeffs::initialise()
{
    if (initialised_)
    {
        return;
    }

    initialised_ = true;
}


void Foam::functionObjects::forceCoeffs::reset()
{
    Cf_.reset();
    Cm_.reset();

    forceCoeff_ == dimensionedVector(dimless, Zero);
    momentCoeff_ == dimensionedVector(dimless, Zero);
}


void Foam::functionObjects::forceCoeffs::calcForceCoeffs()
{
    // Calculate scaling factors
    const scalar pDyn = 0.5*rhoRef_*sqr(magUInf_);
    const dimensionedScalar forceScaling
    (
        dimless/dimForce,
        scalar(1)/(Aref_*pDyn + SMALL)
    );

    // Calculate force coefficients
    forceCoeff_ = forceScaling*force_;

    Cf_.reset
    (
        forceScaling.value()*sumPatchForcesN_,
        forceScaling.value()*sumPatchForcesT_,
        forceScaling.value()*sumInternalForces_
    );
}


void Foam::functionObjects::forceCoeffs::calcMomentCoeffs()
{
    // Calculate scaling factors
    const scalar pDyn = 0.5*rhoRef_*sqr(magUInf_);
    const dimensionedScalar momentScaling
    (
        dimless/(dimForce*dimLength),
        scalar(1)/(Aref_*pDyn*lRef_ + SMALL)
    );

    // Calculate moment coefficients
    momentCoeff_ = momentScaling*moment_;

    Cm_.reset
    (
        momentScaling.value()*sumPatchMomentsN_,
        momentScaling.value()*sumPatchMomentsT_,
        momentScaling.value()*sumInternalMoments_
    );
}


void Foam::functionObjects::forceCoeffs::createIntegratedDataFile()
{
    if (!coeffFilePtr_.valid())
    {
        coeffFilePtr_ = createFile("coefficient");
        writeIntegratedDataFileHeader("Coefficients", coeffFilePtr_());
    }
}


void Foam::functionObjects::forceCoeffs::writeIntegratedDataFileHeader
(
    const word& header,
    OFstream& os
) const
{
    const auto& coordSys = coordSysPtr_();

    writeHeader(os, "Force and moment coefficients");
    writeHeaderValue(os, "dragDir", coordSys.e1());
    writeHeaderValue(os, "sideDir", coordSys.e2());
    writeHeaderValue(os, "liftDir", coordSys.e3());
    writeHeaderValue(os, "rollAxis", coordSys.e1());
    writeHeaderValue(os, "pitchAxis", coordSys.e2());
    writeHeaderValue(os, "yawAxis", coordSys.e3());
    writeHeaderValue(os, "magUInf", magUInf_);
    writeHeaderValue(os, "lRef", lRef_);
    writeHeaderValue(os, "Aref", Aref_);
    writeHeaderValue(os, "CofR", coordSys.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");

    const auto& allCoeffs = coeffsSelection();

    forAllIters(allCoeffs, iter)
    {
        const auto& coeff = iter();

        if (!coeff.active_) continue;

        if (coeff.c1_ == -1)
        {
            writeTabbed(os, coeff.name_);
        }
        else
        {
            writeTabbed(os, coeff.name_);

            if (coeff.splitFrontRear_)
            {
                writeTabbed(os, coeff.frontName());
                writeTabbed(os, coeff.rearName());
            }
        }
    }

    os  << endl;
}


void Foam::functionObjects::forceCoeffs::writeIntegratedDataFile()
{
    OFstream& os = coeffFilePtr_();

    writeCurrentTime(os);

    const vector Cft = Cf_.total();
    const vector Cmt = Cm_.total();

    const auto& allCoeffs = coeffsSelection();

    forAllIters(allCoeffs, iter)
    {
        const auto& coeff = iter();

        if (!coeff.active_) continue;

        if (coeff.c1_ == -1)
        {
            os  << tab << Cmt[coeff.c0_];
        }
        else
        {
            os  << tab << Cft[coeff.c0_];

            if (coeff.splitFrontRear_)
            {
                os  << tab << 0.5*Cft[coeff.c0_] + Cmt[coeff.c1_];
                os  << tab << 0.5*Cft[coeff.c0_] - Cmt[coeff.c1_];
            }
        }
    }

    coeffFilePtr_() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffs::forceCoeffs
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool readFields
)
:
    forces(name, runTime, dict, false),
    Cf_(),
    Cm_(),
    forceCoeff_
    (
        IOobject
        (
            scopedName("forceCoeff"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, Zero)
    ),
    momentCoeff_
    (
        IOobject
        (
            scopedName("momentCoeff"),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, Zero)
    ),
    coeffFilePtr_(),
    magUInf_(Zero),
    lRef_(Zero),
    Aref_(Zero),
    calcFrontRear_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict, "liftDir", "dragDir");
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forceCoeffs::read(const dictionary& dict)
{
    if (!forces::read(dict))
    {
        return false;
    }

    initialised_ = false;

    // Reference scales
    dict.readEntry("magUInf", magUInf_);
    dict.readEntry("lRef", lRef_);
    dict.readEntry("Aref", Aref_);

    // If case is compressible we must read rhoInf
    // (stored in rhoRef_) to calculate the reference dynamic pressure
    // Note: for incompressible, rhoRef_ is already initialised to 1
    if (rhoName_ != "rhoInf")
    {
        dict.readEntry("rhoInf", rhoRef_);
    }

    Info<< "    magUInf: " << magUInf_ << nl
        << "    lRef   : " << lRef_ << nl
        << "    Aref   : " << Aref_ << nl
        << "    rhoInf : " << rhoRef_ << endl;

    if (min(magUInf_, rhoRef_) < SMALL || min(lRef_, Aref_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Non-zero values are required for reference scales"
            << exit(FatalIOError);

        return false;
    }

    if (!dict.found("coefficients"))
    {
        Info<< "Selecting all coefficients" << nl;

        auto& allCoeffs = coeffsSelection();

        forAllIters(allCoeffs, iter)
        {
            iter().active_ = true;
            iter().splitFrontRear_ = true;
            Info<< "- " << iter() << nl;
        }
    }
    else
    {
        wordHashSet coeffs(dict.get<wordHashSet>("coefficients"));
        bool splitFrontRear = dict.get<bool>("splitFrontRear");

        auto& allCoeffs = coeffsSelection();

        Info<< "Selecting coefficients:" << nl;

        forAllIters(coeffs, iter)
        {
            auto& coeff = allCoeffs[iter()];
            coeff.active_ = true;
            coeff.splitFrontRear_ = splitFrontRear;
            Info<< "- " << coeff << nl;
        }
    }

    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;
    }

    Info<< endl;

    return true;
}



bool Foam::functionObjects::forceCoeffs::execute()
{
    forces::calcForcesMoments();

    initialise();

    reset();

    Log << type() << " " << name() << " write:" << nl;

    calcForceCoeffs();

    calcMomentCoeffs();

    auto logValues = [](const word& name, const vector& coeff, Ostream& os)
    {
        os  << "    " << name << ":"
            << tab << "total:" << coeff.x() + coeff.y() + coeff.z()
            << tab << "pressure:" << coeff.x()
            << tab << "viscous:" << coeff.y()
            << tab << "internal:" << coeff.z()
            << nl;
    };

    const vector Cft = Cf_.total();
    const vector Cmt = Cm_.total();

    const auto& allCoeffs = coeffsSelection();

    // Always setting all results
    forAllIters(allCoeffs, iter)
    {
        const auto& coeff = iter();

        if (coeff.c1_ == -1)
        {
            const vector Cmi = Cm_[coeff.c0_];

            if (log && coeff.active_)
            {
                logValues(coeff.name_, Cmi, Info);
            }

            setResult(coeff.name_, Cft[coeff.c0_]);
            setResult(coeff.name_ & "pressure", Cmi.x());
            setResult(coeff.name_ & "viscous", Cmi.y());
            setResult(coeff.name_ & "internal", Cmi.z());
        }
        else
        {
            const vector Cfi = Cf_[coeff.c0_];

            if (log && coeff.active_)
            {
                logValues(coeff.name_, Cfi, Info);

                if (coeff.splitFrontRear_)
                {
                    const vector Cmi = Cm_[coeff.c1_];

                    logValues(coeff.frontName(), 0.5*Cfi + Cmi, Info);
                    logValues(coeff.rearName(), 0.5*Cfi - Cmi, Info);
                }
            }

            setResult(coeff.name_, Cft[coeff.c0_]);
            setResult(coeff.name_ & "pressure", Cfi.x());
            setResult(coeff.name_ & "viscous", Cfi.y());
            setResult(coeff.name_ & "internal", Cfi.z());
            setResult(coeff.frontName(), 0.5*Cft[coeff.c0_] + Cmt[coeff.c1_]);
            setResult(coeff.rearName(), 0.5*Cft[coeff.c0_] - Cmt[coeff.c1_]);
        }
    }

    Log  << endl;

    return true;
}


bool Foam::functionObjects::forceCoeffs::write()
{
    if (writeToFile())
    {
        Log << "    writing force and moment coefficient files." << endl;

        createIntegratedDataFile();

        writeIntegratedDataFile();
    }

    if (writeFields_)
    {
        forceCoeff_.write();
        momentCoeff_.write();
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
