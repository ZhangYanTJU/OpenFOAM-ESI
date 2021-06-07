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


const Foam::Enum
<
    Foam::functionObjects::forceCoeffs::forceCoeffType
>
Foam::functionObjects::forceCoeffs::forceCoeffTypeNames_
({
    { forceCoeffType::CD, "Cd" },
    { forceCoeffType::CS, "Cs" },
    { forceCoeffType::CL, "Cl" }
});


const Foam::Enum
<
    Foam::functionObjects::forceCoeffs::momentCoeffType
>
Foam::functionObjects::forceCoeffs::momentCoeffTypeNames_
({
    { momentCoeffType::CMROLL, "CmRoll" },
    { momentCoeffType::CMPITCH, "CmPitch" },
    { momentCoeffType::CMYAW, "CmYaw" }
});


const Foam::Enum
<
    Foam::functionObjects::forceCoeffs::frontRearForceCoeffType
>
Foam::functionObjects::forceCoeffs::frontRearForceCoeffTypeNames_
({
    { frontRearForceCoeffType::CDF, "Cdf" },
    { frontRearForceCoeffType::CSF, "Csf" },
    { frontRearForceCoeffType::CLF, "Clf" },
    { frontRearForceCoeffType::CDR, "Cdr" },
    { frontRearForceCoeffType::CSR, "Csr" },
    { frontRearForceCoeffType::CLR, "Clr" }
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forceCoeffs::initialise()
{
    if (initialised_)
    {
        return;
    }

    // Collect the identifiers of operand coeffs,
    // which can be different from output coeffs
    idCFs_.clear();
    idCMs_.clear();
    for (const word& coeff : forceCoeffNames_)
    {
        switch (forceCoeffTypeNames_[coeff])
        {
            case CD:
            {
                idCFs_.push_back(0);
                break;
            }
            case CS:
            {
                idCFs_.push_back(1);
                break;
            }
            case CL:
            {
                idCFs_.push_back(2);
                break;
            }
            default:
                break;
        }
    }

    for (const word& coeff : momentCoeffNames_)
    {
        switch (momentCoeffTypeNames_[coeff])
        {
            case CMROLL:
            {
                idCMs_.push_back(0);
                break;
            }
            case CMPITCH:
            {
                idCMs_.push_back(1);
                break;
            }
            case CMYAW:
            {
                idCMs_.push_back(2);
                break;
            }
            default:
                break;
        }
    }

    if (calcFrontRear_)
    {
        for (const word& coeff : frontRearForceCoeffNames_)
        {
            switch (frontRearForceCoeffTypeNames_[coeff])
            {
                case CDF:
                {
                    idCFs_.push_back(0);
                    idCMs_.push_back(0);
                    break;
                }
                case CSF:
                {
                    idCFs_.push_back(1);
                    idCMs_.push_back(2);
                    break;
                }
                case CLF:
                {
                    idCFs_.push_back(2);
                    idCMs_.push_back(1);
                    break;
                }
                case CDR:
                {
                    idCFs_.push_back(0);
                    idCMs_.push_back(0);
                    break;
                }
                case CSR:
                {
                    idCFs_.push_back(1);
                    idCMs_.push_back(2);
                    break;
                }
                case CLR:
                {
                    idCFs_.push_back(2);
                    idCMs_.push_back(1);
                    break;
                }
            }
        }
    }

    // Sort and remove duplicated identifiers of operand coeffs
    std::sort(idCFs_.begin(), idCFs_.end());
    idCFs_.erase
    (
        std::unique(idCFs_.begin(), idCFs_.end()),
        idCFs_.end()
    );
    std::sort(idCMs_.begin(), idCMs_.end());
    idCMs_.erase
    (
        std::unique(idCMs_.begin(), idCMs_.end()),
        idCMs_.end()
    );

    for (label i : idCFs_)
    {
        CFs_[i].setSize(vector::nComponents);
        sumCFs_[i].setSize(vector::nComponents + 1);
    }

    for (label i : idCMs_)
    {
        CMs_[i].setSize(vector::nComponents);
        sumCMs_[i].setSize(vector::nComponents + 1);
    }

    if (calcFrontRear_)
    {
        sumFrontCFs_.setSize(vector::nComponents);
        sumRearCFs_.setSize(vector::nComponents);

        for (direction i = 0; i < vector::nComponents; ++i)
        {
            sumFrontCFs_[i].setSize(vector::nComponents + 1);
            sumRearCFs_[i].setSize(vector::nComponents + 1);
        }
    }

    initialised_ = true;
}


void Foam::functionObjects::forceCoeffs::reset()
{
    CFs_ = Zero;
    CMs_ = Zero;
    sumCFs_ = Zero;
    sumCMs_ = Zero;

    if (calcFrontRear_)
    {
        sumFrontCFs_ = Zero;
        sumRearCFs_ = Zero;
    }

    if (writeFields_)
    {
        auto& fldForceCoeff =
            lookupObjectRef<volVectorField>(fieldName("forceCoeff"));

        fldForceCoeff == dimensionedVector(fldForceCoeff.dimensions(), Zero);

        auto& fldMomentCoeff =
            lookupObjectRef<volVectorField>(fieldName("momentCoeff"));

        fldMomentCoeff == dimensionedVector(fldMomentCoeff.dimensions(), Zero);
    }
}


void Foam::functionObjects::forceCoeffs::calcForceCoeffs()
{
    if (idCFs_.size() == 0)
    {
        return;
    }

    // Calculate scaling factors
    const scalar pDyn = 0.5*rhoRef_*sqr(magUInf_);
    const scalar forceScaling = scalar(1)/(Aref_*pDyn + SMALL);

    // Calculate force coefficients
    for (direction compi = 0; compi < vector::nComponents ; ++compi)
    {
        const Field<vector> localForce(coordSys_.localVector(forces_[compi]));

        for (label coeffi : idCFs_)
        {
            CFs_[coeffi][compi] = forceScaling*(localForce.component(coeffi));
        }
    }

    for (label coeffi : idCFs_)
    {
        // Calculate integrated force coefficients
        sumCFs_[coeffi] = integrateData(CFs_[coeffi]);
    }
}


void Foam::functionObjects::forceCoeffs::calcMomentCoeffs()
{
    if (idCMs_.size() == 0)
    {
        return;
    }

    // Calculate scaling factors
    const scalar pDyn = 0.5*rhoRef_*sqr(magUInf_);
    const scalar momentScaling = scalar(1)/(Aref_*pDyn*lRef_ + SMALL);

    // Calculate moment coefficients
    for (direction compi = 0; compi < vector::nComponents ; ++compi)
    {
        const Field<vector> localMoment(coordSys_.localVector(moments_[compi]));

        for (label coeffi : idCMs_)
        {
            CMs_[coeffi][compi] = momentScaling*(localMoment.component(coeffi));
        }
    }

    for (label coeffi : idCMs_)
    {
        // Calculate integrated moment coefficients
        sumCMs_[coeffi] = integrateData(CMs_[coeffi]);
    }
}


Foam::List<Foam::scalar> Foam::functionObjects::forceCoeffs::integrateData
(
    const List<Field<scalar>>& coeff
) const
{
    scalar pressure = sum(coeff[0]);
    scalar viscous = sum(coeff[1]);
    scalar porous = porosity_ ? sum(coeff[2]) : 0;
    scalar total = pressure + viscous + porous;

    return List<scalar>({total, pressure, viscous, porous});
}


bool Foam::functionObjects::forceCoeffs::createDataFile()
{
    if (!coeffFilePtr_.valid())
    {
        coeffFilePtr_ = createFile("coefficient");
        writeDataFileHeader("Coefficients", coeffFilePtr_());
    }

    return true;
}


bool Foam::functionObjects::forceCoeffs::writeDataFileHeader
(
    const word& header,
    OFstream& os
) const
{
    writeHeader(os, "Force and moment coefficients");
    writeHeaderValue(os, "dragDir", coordSys_.e1());
    writeHeaderValue(os, "sideDir", coordSys_.e2());
    writeHeaderValue(os, "liftDir", coordSys_.e3());
    writeHeaderValue(os, "rollAxis", coordSys_.e1());
    writeHeaderValue(os, "pitchAxis", coordSys_.e2());
    writeHeaderValue(os, "yawAxis", coordSys_.e3());
    writeHeaderValue(os, "magUInf", magUInf_);
    writeHeaderValue(os, "lRef", lRef_);
    writeHeaderValue(os, "Aref", Aref_);
    writeHeaderValue(os, "CofR", coordSys_.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");

    for (const word& coeff : forceCoeffNames_)
    {
        writeTabbed(os, coeff);
    }

    for (const word& coeff : momentCoeffNames_)
    {
        writeTabbed(os, coeff);
    }

    if (calcFrontRear_)
    {
        for (const word& coeff : frontRearForceCoeffNames_)
        {
            writeTabbed(os, coeff);
        }
    }

    os  << endl;

    return true;
}


bool Foam::functionObjects::forceCoeffs::writeDataFile()
{
    writeCurrentTime(coeffFilePtr_());

    for (const word& coeff : forceCoeffNames_)
    {
        switch (forceCoeffTypeNames_[coeff])
        {
            case CD:
            {
                coeffFilePtr_() << sumCFs_[0][0] << tab;
                break;
            }
            case CS:
            {
                coeffFilePtr_() << sumCFs_[1][0] << tab;
                break;
            }
            case CL:
            {
                coeffFilePtr_() << sumCFs_[2][0] << tab;
                break;
            }
            default:
                break;
        }
    }

    for (const word& coeff : momentCoeffNames_)
    {
        switch (momentCoeffTypeNames_[coeff])
        {
            case CMROLL:
            {
                coeffFilePtr_() << sumCMs_[0][0] << tab;
                break;
            }
            case CMPITCH:
            {
                coeffFilePtr_() << sumCMs_[1][0] << tab;
                break;
            }
            case CMYAW:
            {
                coeffFilePtr_() << sumCMs_[2][0] << tab;
                break;
            }
            default:
                break;
        }
    }

    if (calcFrontRear_)
    {
        for (const word& coeff : frontRearForceCoeffNames_)
        {
            switch (frontRearForceCoeffTypeNames_[coeff])
            {
                case CDF:
                {
                    coeffFilePtr_() << sumFrontCFs_[0][0] << tab;
                    break;
                }
                case CSF:
                {
                    coeffFilePtr_() << sumFrontCFs_[1][0] << tab;
                    break;
                }
                case CLF:
                {
                    coeffFilePtr_() << sumFrontCFs_[2][0] << tab;
                    break;
                }
                case CDR:
                {
                    coeffFilePtr_() << sumRearCFs_[0][0] << tab;
                    break;
                }
                case CSR:
                {
                    coeffFilePtr_() << sumRearCFs_[1][0] << tab;
                    break;
                }
                case CLR:
                {
                    coeffFilePtr_() << sumRearCFs_[2][0] << tab;
                    break;
                }
            }
        }
    }

    coeffFilePtr_() <<  endl;

    return true;
}


void Foam::functionObjects::forceCoeffs::logData
(
    const word& title,
    const List<scalar>& coeff
) const
{
    if (!log)
    {
        return;
    }

    Log << "        " << title
        << "       : "  << coeff[0] << tab
        << "pressure: " << coeff[1] << tab
        << "viscous: "  << coeff[2];

    if (porosity_)
    {
        Log << tab << "porous: " << coeff[3];
    }

    Log << endl;
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
    CFs_(vector::nComponents),
    CMs_(vector::nComponents),
    sumCFs_(vector::nComponents),
    sumCMs_(vector::nComponents),
    sumFrontCFs_(Zero),
    sumRearCFs_(Zero),
    forceCoeffNames_(),
    momentCoeffNames_(),
    frontRearForceCoeffNames_(),
    idCFs_(),
    idCMs_(),
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

    forceCoeffNames_.insert
    (
        dict.getOrDefault<wordList>
        (
            "forceCoeffs",
            wordList({"Cd", "Cs", "Cl"})
        )
    );

    if (forceCoeffNames_.size())
    {
        Info<< "    Operand force coefficients:" << nl;
        for (const word& forceCoeff : forceCoeffNames_)
        {
            Info<< "        " << forceCoeff << nl;
        }
    }

    momentCoeffNames_.insert
    (
        dict.getOrDefault<wordList>
        (
            "momentCoeffs",
            wordList({"CmRoll", "CmPitch", "CmYaw"})
        )
    );

    if (momentCoeffNames_.size())
    {
        Info<< "    Operand moment coefficients:" << nl;
        for (const word& momentCoeff : momentCoeffNames_)
        {
            Info<< "        " << momentCoeff << nl;
        }
    }

    if (dict.found("frontRearForceCoeffs"))
    {
        frontRearForceCoeffNames_.insert
        (
            dict.get<wordList>("frontRearForceCoeffs")
        );
    }

    if (frontRearForceCoeffNames_.size())
    {
        Info<< "    Operand front- and rear-axle force coefficients:" << nl;
        for (const word& frontRearForceCoeff : frontRearForceCoeffNames_)
        {
            Info<< "        " << frontRearForceCoeff << nl;
        }

        calcFrontRear_ = true;
    }

    if
    (
        (
            forceCoeffNames_.size()
          + momentCoeffNames_.size()
          + frontRearForceCoeffNames_.size()
        )
        <= 0
    )
    {
        FatalIOErrorInFunction(dict)
            << "No coefficient is requested" << nl
            << "Please request at least one force/moment coefficient"
            << exit(FatalIOError);

        return false;
    }

    if (writeFields_)
    {
        auto* fldForceCoeffPtr
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("forceCoeff"),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector(dimless, Zero)
            )
        );

        mesh_.objectRegistry::store(fldForceCoeffPtr);

        auto* fldMomentCoeffPtr
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("momentCoeff"),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector(dimless, Zero)
            )
        );

        mesh_.objectRegistry::store(fldMomentCoeffPtr);
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

    // Calculate force coefficients
    for (const word& coeff : forceCoeffNames_)
    {
        switch (forceCoeffTypeNames_[coeff])
        {
            case CD:
            {
                setResult(coeff, sumCFs_[0][0]);
                logData(coeff, sumCFs_[0]);
                break;
            }
            case CS:
            {
                setResult(coeff, sumCFs_[1][0]);
                logData(coeff, sumCFs_[1]);
                break;
            }
            case CL:
            {
                setResult(coeff, sumCFs_[2][0]);
                logData(coeff, sumCFs_[2]);
                break;
            }
            default:
                break;
        }
    }

    // Calculate moment coefficients
    for (const word& coeff : momentCoeffNames_)
    {
        switch (momentCoeffTypeNames_[coeff])
        {
            case CMROLL:
            {
                setResult(coeff, sumCMs_[0][0]);
                logData(coeff, sumCMs_[0]);
                break;
            }
            case CMPITCH:
            {
                setResult(coeff, sumCMs_[1][0]);
                logData(coeff, sumCMs_[1]);
                break;
            }
            case CMYAW:
            {
                setResult(coeff, sumCMs_[2][0]);
                logData(coeff, sumCMs_[2]);
                break;
            }
            default:
                break;
        }
    }

    if (calcFrontRear_)
    {
        // Calculate integrated front- and rear-axle force coefficients
        for (const word& coeff : frontRearForceCoeffNames_)
        {
            switch (frontRearForceCoeffTypeNames_[coeff])
            {
                case CDF:
                {
                    sumFrontCFs_[0] = integrateData(0.5*CFs_[0] + CMs_[0]);
                    setResult(coeff, sumFrontCFs_[0][0]);
                    logData(coeff, sumFrontCFs_[0]);
                    break;
                }
                case CSF:
                {
                    sumFrontCFs_[1] = integrateData(0.5*CFs_[1] + CMs_[2]);
                    setResult(coeff, sumFrontCFs_[1][0]);
                    logData(coeff, sumFrontCFs_[1]);
                    break;
                }
                case CLF:
                {
                    sumFrontCFs_[2] = integrateData(0.5*CFs_[2] + CMs_[1]);
                    setResult(coeff, sumFrontCFs_[2][0]);
                    logData(coeff, sumFrontCFs_[2]);
                    break;
                }
                case CDR:
                {
                    sumRearCFs_[0] = integrateData(0.5*CFs_[0] - CMs_[0]);
                    setResult(coeff, sumRearCFs_[0][0]);
                    logData(coeff, sumRearCFs_[0]);
                    break;
                }
                case CSR:
                {
                    sumRearCFs_[1] = integrateData(0.5*CFs_[1] - CMs_[2]);
                    setResult(coeff, sumRearCFs_[1][0]);
                    logData(coeff, sumRearCFs_[1]);
                    break;
                }
                case CLR:
                {
                    sumRearCFs_[2] = integrateData(0.5*CFs_[2] - CMs_[1]);
                    setResult(coeff, sumRearCFs_[2][0]);
                    logData(coeff, sumRearCFs_[2]);
                    break;
                }
            }
        }
    }


    if (writeFields_)
    {
        const auto& fldForce =
            lookupObject<volVectorField>(fieldName("force"));

        const auto& fldMoment =
            lookupObject<volVectorField>(fieldName("moment"));

        auto& fldForceCoeff =
            lookupObjectRef<volVectorField>(fieldName("forceCoeff"));

        auto& fldMomentCoeff =
            lookupObjectRef<volVectorField>(fieldName("momentCoeff"));

        const scalar pDyn = 0.5*rhoRef_*sqr(magUInf_);
        const dimensionedScalar f0(dimForce, Aref_*pDyn);
        const dimensionedScalar m0(dimForce*dimLength, Aref_*lRef_*pDyn);

        fldForceCoeff == fldForce/f0;
        fldMomentCoeff == fldMoment/m0;
    }

    return true;
}


bool Foam::functionObjects::forceCoeffs::write()
{
    if (writeToFile() && Pstream::master())
    {
        Log << "    writing force and moment coefficient files." << endl;

        createDataFile();

        writeDataFile();
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
