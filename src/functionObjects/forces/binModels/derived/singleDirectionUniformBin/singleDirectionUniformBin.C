/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "singleDirectionUniformBin.H"
#include "porosityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{
    defineTypeNameAndDebug(singleDirectionUniformBin, 0);
    addToRunTimeSelectionTable(binModel, singleDirectionUniformBin, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binModels::singleDirectionUniformBin::singleDirectionUniformBin
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    binModel(name, dict, mesh),
    writeFile(mesh, name),
    cumulative_(true),
    binDx_(0),
    binMin_(GREAT),
    binMax_(GREAT),
    binDir_(Zero)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::binModels::singleDirectionUniformBin::read(const dictionary& dict)
{
    if (!binModel::read(dict))
    {
        return false;
    }

    if (dict.found("binData"))
    {
        Info<< "    Activating a set of single-direction bins" << endl;

        const dictionary& binDict(dict.subDict("binData"));

        nBin_ = binDict.get<label>("nBin");

        if (nBin_ <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "Number of bins must be greater than zero" << nl
                << "    nBin = " << nBin_ << nl
                << exit(FatalIOError);
        }
        else
        {
            Info<< "    Employing " << nBin_ << " bins" << endl;
            if (binDict.readIfPresent("min", binMin_))
            {
                Info<< "    - min        : " << binMin_ << endl;
            }
            if (binDict.readIfPresent("max", binMax_))
            {
                Info<< "    - max        : " << binMax_ << endl;
            }

            cumulative_ = binDict.get<bool>("cumulative");
            Info<< "    - cumulative    : " << cumulative_ << endl;

            binDir_ = binDict.get<vector>("direction");
            binDir_.normalise();

            if (mag(binDir_) == 0)
            {
                FatalIOErrorInFunction(dict)
                    << "Input direction should not be zero valued" << nl
                    << "    direction = " << binDir_ << nl
                    << exit(FatalIOError);
            }

            Info<< "    - direction     : " << binDir_ << endl;
        }
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "No bin properties are given by 'binData' dictionary" << nl
            << exit(FatalIOError);

        return false;
    }

    return true;
}


void Foam::binModels::singleDirectionUniformBin::initialise()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Determine extents of patches in a given direction
    scalar geomMin = GREAT;
    scalar geomMax = -GREAT;
    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = pbm[patchi];
        const scalarField d(pp.faceCentres() & binDir_);
        geomMin = min(min(d), geomMin);
        geomMax = max(max(d), geomMax);
    }

    if (porosity_)
    {
        const HashTable<const porosityModel*> models =
            mesh_.lookupClass<porosityModel>();

        const scalarField dd(mesh_.C() & binDir_);

        forAllConstIters(models, iter)
        {
            const porosityModel& pm = *iter();
            const labelList& cellZoneIDs = pm.cellZoneIDs();

            for (const label zonei : cellZoneIDs)
            {
                const cellZone& cZone = mesh_.cellZones()[zonei];
                const scalarField d(dd, cZone);
                geomMin = min(min(d), geomMin);
                geomMax = max(max(d), geomMax);
            }
        }
    }

    reduce(geomMin, minOp<scalar>());
    reduce(geomMax, maxOp<scalar>());

    // Slightly boost max so that region of interest is fully within bounds
    geomMax = 1.0001*(geomMax - geomMin) + geomMin;

    // Use geometry limits if not specified by the user
    if (binMin_ == GREAT) binMin_ = geomMin;
    if (binMax_ == GREAT) binMax_ = geomMax;

    binDx_ = (binMax_ - binMin_)/scalar(nBin_);

    if (binDx_ <= 0)
    {
        FatalErrorInFunction
            << "Max bound must be greater than min bound" << nl
            << "    d           = " << binDx_ << nl
            << "    min         = " << binMin_ << nl
            << "    max         = " << binMax_ << nl
            << exit(FatalError);
    }
}


void Foam::binModels::singleDirectionUniformBin::applyBins
(
    List<Field<vector>>& force,
    List<Field<vector>>& moment,
    const vectorField& d,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT
)
{
    const scalarField dd(d & binDir_);

    forAll(dd, i)
    {
        // Avoid faces outside of the bin
        bool faceInside = true;
        if (dd[i] < binMin_ || dd[i] > binMax_)
        {
            faceInside = false;
        }

        if (faceInside)
        {
            // Find the bin division corresponding to the face
            const label bini =
                min(max(floor((dd[i] - binMin_)/binDx_), 0), nBin_ - 1);

            force[0][bini] += fN[i];
            force[1][bini] += fT[i];
            moment[0][bini] += Md[i]^fN[i];
            moment[1][bini] += Md[i]^fT[i];
        }
    }
}


void Foam::binModels::singleDirectionUniformBin::applyBins
(
    List<Field<vector>>& force,
    List<Field<vector>>& moment,
    const vectorField& d,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    const scalarField dd(d & binDir_);

    forAll(dd, i)
    {
        // Avoid faces outside of the bin
        bool faceInside = true;
        if (dd[i] < binMin_ || dd[i] > binMax_)
        {
            faceInside = false;
        }

        if (faceInside)
        {
            // Find the bin division corresponding to the face
            const label bini =
                min(max(floor((dd[i] - binMin_)/binDx_), 0), nBin_ - 1);

            force[0][bini] += fN[i];
            force[1][bini] += fT[i];
            force[2][bini] += fP[i];
            moment[0][bini] += Md[i]^fN[i];
            moment[1][bini] += Md[i]^fT[i];
            moment[2][bini] += Md[i]^fP[i];
        }
    }
}


void Foam::binModels::singleDirectionUniformBin::createBinnedDataFiles()
{
    if (!forceBinFilePtr_.valid())
    {
        forceBinFilePtr_ = createFile("forceBin");
        writeBinnedDataFileHeader("Force", forceBinFilePtr_());
    }

    if (!momentBinFilePtr_.valid())
    {
        momentBinFilePtr_ = createFile("momentBin");
        writeBinnedDataFileHeader("Moment", momentBinFilePtr_());
    }
}


void Foam::binModels::singleDirectionUniformBin::createBinnedDataFile
(
    const word& title,
    const label id
)
{
    if (!coeffBinFilePtrs_.set(id))
    {
        coeffBinFilePtrs_.set(id, createFile(title + "Bin"));
        OFstream& os = coeffBinFilePtrs_[id];
        writeBinnedDataFileHeader(title, os, true);
    }
}


void Foam::binModels::singleDirectionUniformBin::writeBinnedDataFileHeader
(
    const word& header,
    OFstream& os,
    const bool coeffStyle
) const
{
    writeHeader(os, header + " bins");
    writeHeaderValue(os, "bins", nBin_);
    writeHeaderValue(os, "start", binMin_);
    writeHeaderValue(os, "end", binMax_);
    writeHeaderValue(os, "delta", binDx_);
    writeHeaderValue(os, "direction", binDir_);

    // Compute and print bin end points in the binning direction
    vectorField binPoints(nBin_);
    writeCommented(os, "x co-ords  :");
    forAll(binPoints, pointi)
    {
        binPoints[pointi] = (binMin_ + (pointi + 1)*binDx_)*binDir_;
        os  << tab << binPoints[pointi].x();
    }
    os  << nl;

    writeCommented(os, "y co-ords  :");
    forAll(binPoints, pointi)
    {
        os  << tab << binPoints[pointi].y();
    }
    os  << nl;

    writeCommented(os, "z co-ords  :");
    forAll(binPoints, pointi)
    {
        os  << tab << binPoints[pointi].z();
    }
    os  << nl;

    writeHeader(os, "");
    writeCommented(os, "Time");

    if (coeffStyle)
    {
        for (label i = 0; i < nBin_; ++i)
        {
            const word in(Foam::name(i) + ':');
            writeTabbed(os, in + "total");
            writeTabbed(os, in + "pressure");
            writeTabbed(os, in + "viscous");

            if (porosity_)
            {
                writeTabbed(os, in + "porous");
            }
        }
    }
    else
    {
        for (label i = 0; i < nBin_; ++i)
        {
            const word in(Foam::name(i) + ':');
            writeTabbed(os, in + "(total_x total_y total_z)");
            writeTabbed(os, in + "(pressure_x pressure_y pressure_z)");
            writeTabbed(os, in + "(viscous_x viscous_y viscous_z)");

            if (porosity_)
            {
                writeTabbed(os, in + "(porous_x porous_y porous_z)");
            }
        }
    }

    os  << endl;
}


void Foam::binModels::singleDirectionUniformBin::writeBinnedDataFiles
(
    const List<Field<vector>>& force,
    const List<Field<vector>>& moment,
    const coordSystem::cartesian& coordSys
)
{
    List<Field<vector>> lf(vector::nComponents);
    List<Field<vector>> lm(vector::nComponents);

    for (direction i = 0; i < vector::nComponents; ++i)
    {
        lf[i] = coordSys.localVector(force[i]);
        lm[i] = coordSys.localVector(moment[i]);
    }

    writeBinnedDataFile(lf, forceBinFilePtr_);
    writeBinnedDataFile(lm, momentBinFilePtr_);
}


void Foam::binModels::singleDirectionUniformBin::writeBinnedDataFile
(
    List<Field<vector>>& fm,
    autoPtr<OFstream>& osPtr
) const
{
    if (cumulative_)
    {
        for (direction i = 0; i < vector::nComponents; ++i)
        {
            for (label bini = 1; bini < nBin_; ++bini)
            {
                fm[i][bini] += fm[i][bini-1];
            }
        }
    }

    OFstream& os = osPtr();

    writeCurrentTime(os);

    for (label bini = 0; bini < nBin_; ++bini)
    {
        const vector total(fm[0][bini] + fm[1][bini] + fm[2][bini]);

        os  << tab << total
            << tab << fm[0][bini]
            << tab << fm[1][bini];

        if (porosity_)
        {
            os  << tab << fm[2][bini];
        }
    }

    os  << nl;
}


void Foam::binModels::singleDirectionUniformBin::writeBinnedDataFile
(
    List<Field<scalar>> coeff,
    const label id
)
{
    if (cumulative_)
    {
        for (direction i = 0; i < vector::nComponents; ++i)
        {
            for (label bini = 1; bini < nBin_; ++bini)
            {
                coeff[i][bini] += coeff[i][bini-1];
            }
        }
    }

    OFstream& os = coeffBinFilePtrs_[id];

    writeCurrentTime(os);

    for (label bini = 0; bini < nBin_; ++bini)
    {
        scalar total = coeff[0][bini] + coeff[1][bini];

        if (porosity_)
        {
            total += coeff[2][bini];
        }

        os  << tab << total
            << tab << coeff[0][bini]
            << tab << coeff[1][bini];

        if (porosity_)
        {
            os  << tab << coeff[2][bini];
        }
    }

    os  << endl;
}


// ************************************************************************* //
