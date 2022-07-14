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

#include "compressedIncompressibleVars.H"
#include "createZeroField.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(compressedIncompressibleVars, 0);
defineRunTimeSelectionTable(compressedIncompressibleVars, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void compressedIncompressibleVars::writeLog
(
    const word& name
) const
{
    const scalar& initialSize = storageMetrics_[0];
    const scalar& uncompressedRoughSize = storageMetrics_[1];
    const scalar& uncompressedSize = storageMetrics_[2];
    const scalar& compressedSize = storageMetrics_[3];
    if (!initialSize || !uncompressedSize || !uncompressedRoughSize)
    {
        FatalErrorInFunction
            << "initialSize, uncompressedRoughSize or uncompressedSize is Zero"
            << endl
            << exit(FatalError);
    }
    if (!compressedSize)
    {
        FatalErrorInFunction
            << "compressedSize is Zero" << endl
            << exit(FatalError);
    }
    unsigned int width = 6;//IOstream::defaultPrecision();
    Info<< setprecision(width);
    Info<< name << tab << "pure algorithm's CR = "
        << uncompressedSize/compressedSize << endl;
    if (debug)
    {
        Info<< name << ": Initial Size             = "
            << initialSize*1.e-6 << " Mb" << nl
            << name << ": Initial Size saved       = "
            << uncompressedRoughSize*1.e-6 << " Mb" << nl
            << name << ": Uncompressed Size        = "
            << uncompressedSize*1.e-6 << " Mb" << nl
            << name << ": Compressed Size          = "
            << compressedSize*1.e-6 << " Mb" << endl;
        if (name == "All")
        {
            Info<< name
                << ": Data Increase due to data manipulation before compression = "
                << (uncompressedSize-uncompressedRoughSize)/uncompressedRoughSize*100. << "%"
                << endl;
        }
        Info<< name << ": (%) Data Reduction       = "
            << (uncompressedSize-compressedSize)/uncompressedSize*100. << "%"
            << endl;
    }
    Info<< setprecision(IOstream::defaultPrecision());
    if (name == "All")
    {
        Info << endl;
    }
}


void compressedIncompressibleVars::write(const label i)
{
    if (Pstream::master())
    {
        if (storageFilePtr_.empty())
        {
            setStorageFilesPtr();
        }
        const scalar initialSize = storageMetrics_[0];
        const scalar uncompressedRoughSize = storageMetrics_[1];
        const scalar uncompressedSize = storageMetrics_[2];
        const scalar compressedSize = storageMetrics_[3];

        unsigned int width = IOstream::defaultPrecision() + 6;
        storageFilePtr_[i]
            << setprecision(6) << setw(8) << mesh_.time().value() << tab
            << setprecision(IOstream::defaultPrecision())
            << setw(width) << uncompressedSize/compressedSize << tab
            << setw(width) << (uncompressedSize-compressedSize)/uncompressedSize*100. << tab
            << setw(width) << initialSize*1.e-6 << tab
            << setw(width) << uncompressedRoughSize*1.e-6 << tab
            << setw(width) << compressedSize*1.e-6 << endl;
        storageFilePtr_.clear();
    }
}


void compressedIncompressibleVars::setFields()
{
    incompressibleVars& v = incoVars_;
    p_.reset(compressedVolScalarField::New(v.pInst(), storageParams_, 0, k_));
    phi_.reset
    (
        compressedSurfaceScalarField::New(v.phiInst(), storageParams_, 1, k_)
    );
    U_.reset(compressedVolVectorField::New(v.UInst(), storageParams_, 2, k_));
    incompressible::RASModelVariables& rasVars = v.RASModelVariables()();
    if (rasVars.hasTMVar1())
    {
        RASModelVars_.append
        (
            compressedVolScalarField::New
                (rasVars.TMVar1Inst(), storageParams_, 3, k_)
        );
    }
    if (rasVars.hasTMVar2())
    {
        RASModelVars_.append
        (
            compressedVolScalarField::New
                (rasVars.TMVar2Inst(), storageParams_, 4, k_)
        );
    }
}


void compressedIncompressibleVars::makeFolder()
{
    if (Pstream::master())
    {
        storageFolder_ = mesh_.time().globalPath()/"optimisation"/"storage";
        if (!isDir(storageFolder_))
        {
            mkDir(storageFolder_);
            WarningInFunction
                << "The folder storing the compression metrics should have "
                << "been already set."
                << nl
                << "Possible error in: compressedFullStorage::makeFolder()."
                << "Setting folder anyway." << nl
                << endl;
        }
    }
}


void compressedIncompressibleVars::setStorageFilesPtr()
{
    storageFilePtr_.setSize(7);
    label Ic = -1;
    names_.setSize(7, "undefined");
    names_[++Ic] = "AllFields" + solverName();
    names_[++Ic] = incoVars_.pInst().name();
    names_[++Ic] = incoVars_.phiInst().name();
    names_[++Ic] = incoVars_.UInst().name();
    incompressible::RASModelVariables& rasVars =
        incoVars_.RASModelVariables()();
    forAll(RASModelVars_, iPtr)
    {
        if (iPtr == 0)
        {
            names_[++Ic] = rasVars.TMVar1Inst().name();
        }
        else if (iPtr == 1)
        {
            names_[++Ic] = rasVars.TMVar2Inst().name();
        }
    }
    for (label i = 0; i <= Ic; i++)
    {
        fileName filePath = storageFolder_/names_[i];
        if (!isFile(filePath))
        {
            storageFilePtr_.set
            (
                i, new OFstream(filePath)
            );
            unsigned int width = IOstream::defaultPrecision() + 6;
            storageFilePtr_[i]
                << setw(8) << "# Time" << tab
                << setw(width) << "algorithm's CR" << tab
                << setw(width) << "(%) Reduction" << tab
                << setw(width) << "Initial Size (Mb)" << tab
                << setw(width) << "Uncompressed Size (Mb)" << tab
                << setw(width) << "Compressed Size (Mb)" << endl;
        }
        else
        {
            storageFilePtr_.set
            (
                i,
                new OFstream
                (
                    filePath,
                    IOstream::ASCII,
                    IOstream::currentVersion,
                    IOstream::UNCOMPRESSED,
                    true
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

compressedIncompressibleVars::compressedIncompressibleVars
(
    incompressibleVars& vs,
    storageParameters& storageParams
)
:
    variablesSet
    (
        vs.mesh(),
        vs.solverControlReference().solverDict()
    ),
    incoVars_(vs),
    solverControl_(vs.solverControlReference()),
    storageParams_(storageParams),
    k_(0),
    kCopy_(k_),
    p_(nullptr),
    phi_(nullptr),
    U_(nullptr),
    RASModelVars_(),
    timeIndex_(mesh_.time().timeIndex()),
    timeValue_(mesh_.time().value()),
    storageMetrics_(),
    names_(7, "undefined"),
    storageFilePtr_()
{
    setFields();
    makeFolder();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<compressedIncompressibleVars> compressedIncompressibleVars::New
(
    incompressibleVars& vs,
    storageParameters& storageParams
)
{
    const word type = storageParams.algorithm() + "IncompressibleVars";
    DebugInfo
        << "compressedIncompressibleVars type for the primal fields: " << type
        << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(type);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown compressedIncompressibleVars type " << type << nl << nl
            << "Valid compressedIncompressibleVars types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<compressedIncompressibleVars>
        (cstrIter()(vs, storageParams));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void compressedIncompressibleVars::calculateStorageMetrics()
{
    storageMetrics_.setSize(p_().storageMetrics().size());
    label i = 0;
    if (!storageParams_.timing())
    {
        calculateAndWrite(p_(), incoVars_.pInst().name(), i);
        calculateAndWrite(phi_(), incoVars_.phiInst().name(), i);
        calculateAndWrite(U_(), incoVars_.UInst().name(), i);
        incompressible::RASModelVariables& rasVars =
            incoVars_.RASModelVariables()();
        if (rasVars.hasTMVar1())
        {
            calculateAndWrite(RASModelVars_[0], rasVars.TMVar1Inst().name(), i);
        }
        if (rasVars.hasTMVar2())
        {
            calculateAndWrite(RASModelVars_[1], rasVars.TMVar2Inst().name(), i);
        }
    }
    storageMetrics_ =
        p_().storageMetrics()
      + phi_().storageMetrics()
      + U_().storageMetrics();

    forAll(RASModelVars_, iPtr)
    {
        storageMetrics_ =
            storageMetrics_ + RASModelVars_[iPtr].storageMetrics();
    }
    if (!storageParams_.timing())
    {
        write(0);
        writeLog("All");
    }
    storageFilePtr_.clear();
    names_.clear();
}


void compressedIncompressibleVars::compress()
{
    NotImplemented
}


void compressedIncompressibleVars::decompress
(
    incompressibleVars& vars
)
{
    NotImplemented
}


void compressedIncompressibleVars::decompress()
{
    this->decompress(incoVars_);
}


const scalarList& compressedIncompressibleVars::storageMetrics() const
{
    return storageMetrics_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
