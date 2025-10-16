/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

#include "radiometerProbes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(radiometerProbes, 0);
    addToRunTimeSelectionTable(functionObject, radiometerProbes, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::radiometerProbes::writeFileHeader(Ostream& os)
{
    const pointField& locs = probeLocations();

    writeCommented(os, "Probe,Location,Normal");

    os  << nl;
    for (label i = 0; i < szProbes_; ++i)
    {
        const vector& loc = locs[i];
        const vector& n = n_[i];

        os  << '#' << ' ' << i
            << ',' << loc.x() << ',' << loc.y() << ',' << loc.z()
            << ',' << n.x()   << ',' << n.y()   << ',' << n.z()
            << nl;
    }

    os  << "# Time";
    for (int i = 0; i < szProbes_; ++i)
    {
        os  << ',' << i;
    }
    os  << endl;

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::radiometerProbes::radiometerProbes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    internalProber(mesh_, dict),
    writeFile(mesh_, name, typeName, dict),
    dom_(mesh_.lookupObject<radiation::fvDOM>("radiationProperties")),
    firstIter_(true)
{
    read(dict);
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::radiometerProbes::read(const dictionary& dict)
{
    if
    (
        !(
            regionFunctionObject::read(dict)
         && internalProber::read(dict)
         && writeFile::read(dict)
        )
    )
    {
        return false;
    }


    // Skip if the radiation model is inactive
    if (!dom_.radiation())
    {
        WarningInFunction
            << "The radiation model is inactive."
            << "Skipping the function object " << type() << ' ' << name()
            << endl;
        return false;
    }

    Log << type() << ':' << name() << ": read" << nl << nl;


    // Probe locations are read by 'internalProber'
    szProbes_ = this->size();

    // If/when fvDOM is updated, the 'read' func is assumed to be executed to
    // update the fvDOM properties
    nRay_ = dom_.nRay();

    if (!szProbes_ || !nRay_)
    {
        FatalIOErrorInFunction(dict)
            << "size(probe locations): " << szProbes_ << nl
            << "size(rays): " << nRay_ << nl
            << "The input size of probe locations and rays cannot be zero."
            << exit(FatalIOError);
    }

    // Read and check size consistency of probe normals with probe locations
    dict.readEntry("probeNormals", n_);
    if (n_.size() != szProbes_)
    {
        FatalIOErrorInFunction(dict)
            << "size(probe locations): " << szProbes_ << nl
            << "size(probe normals): " << n_.size() << nl
            << "The input size of probe locations and normals must match."
            << exit(FatalIOError);
    }
    n_.normalise();


    // Pre-compute and cache inner product of 'n_' and 'dAve_', and 'C_'
    // This simplification of calculation is valid only if I is non-negative
    n_dAve_.resize(nRay_);
    C_.resize(nRay_);

    for (label rayi = 0; rayi < nRay_; ++rayi)
    {
        const vector& dAvei = dom_.IRay(rayi).dAve();

        scalarList& n_dAveRay = n_dAve_[rayi];
        boolList& Cray = C_[rayi];

        n_dAveRay.resize(szProbes_, Zero);
        Cray.resize(szProbes_, false);

        for (label pi = 0; pi < szProbes_; ++pi)
        {
            n_dAveRay[pi] = n_[pi] & dAvei;

            if (n_dAveRay[pi] < 0)  // ray entering the probe
            {
                Cray[pi] = true;
            }
        }
    }

    qin_.resize(szProbes_);


    if (writeFile::canResetFile())
    {
        writeFile::resetFile(typeName);
    }

    if (writeFile::canWriteHeader())
    {
        writeFileHeader(file());
    }

    return true;
}


bool Foam::functionObjects::radiometerProbes::execute()
{
    // Skip if there is no probe to sample, or the radiation model is inactive
    if (!szProbes_ || !dom_.radiation() || !shouldCalcThisStep())
    {
        return false;
    }

    Log << type() << ' ' << name() << ": execute" << nl << nl;

    qin_ = Zero;  // resized in 'read'

    for (label rayi = 0; rayi < nRay_; ++rayi)
    {
        // Radiative intensity field for this ray
        const volScalarField& I = dom_.IRay(rayi).I();

        // Sample radiative intensity ray at probe locations
        tmp<scalarField> tIp  = internalProber::sample(I);
        const scalarField& Ip = tIp.cref();

        const scalarList& n_dAveRay = n_dAve_[rayi];
        const boolList& Cray = C_[rayi];

        // Add incident radiative heat flux per probe location for each ray
        for (label pi = 0; pi < szProbes_; ++pi)
        {
            if (Cray[pi])
            {
                qin_[pi] += Ip[pi]*n_dAveRay[pi];
            }
        }
    }

    return true;
}


bool Foam::functionObjects::radiometerProbes::write()
{
    // Skip if there is no probe to sample, or the radiation model is inactive
    if (!szProbes_ || !dom_.radiation() || !shouldCalcThisStep())
    {
        return false;
    }

    Log << type() << ' ' << name() << ": write" << nl << nl;

    if (UPstream::master())
    {
        Ostream& os = file();

        os  << mesh_.time().timeOutputValue();
        for (label pi = 0; pi < szProbes_; ++pi)
        {
            os  << ',' << qin_[pi];
        }
        os  << endl;
    }

    firstIter_ = false;

    return true;
}


// ************************************************************************* //
