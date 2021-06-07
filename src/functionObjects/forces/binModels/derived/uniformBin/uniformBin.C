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

#include "uniformBin.H"
#include "porosityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{
    defineTypeNameAndDebug(uniformBin, 0);
    addToRunTimeSelectionTable(binModel, uniformBin, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binModels::uniformBin::uniformBin
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    binModel(name, dict, mesh),
    writeFile(mesh, name),
    cumulative_(false),
    nBins_(Zero),
    binW_(Zero),
    binMinMax_
    (
        vector2D(GREAT, GREAT),
        vector2D(GREAT, GREAT),
        vector2D(GREAT, GREAT)
    ),
    csysPtr_()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::binModels::uniformBin::read(const dictionary& dict)
{
    if (!binModel::read(dict))
    {
        return false;
    }

    if (dict.found("binData"))
    {
        Info<< "    Activating a set of uniform bins" << endl;

        const dictionary& binDict(dict.subDict("binData"));

        if (binDict.found(coordinateSystem::typeName_()))
        {
            csysPtr_ =
                coordinateSystem::New
                (
                    mesh_,
                    binDict,
                    coordinateSystem::typeName_()
                );

            Info<< "    - type          : " << csysPtr_->name() << nl
                << "    - origin        : " << csysPtr_->origin() << nl
                << "    - e3            : " << csysPtr_->e3() << nl
                << "    - e1            : " << csysPtr_->e1() << endl;

            nBins_ = binDict.get<Vector<label>>("nBin");
        }
        else
        {
            FatalIOErrorInFunction(binDict)
                << "Missing orientation entries in 'binData' dictionary" << nl
                << exit(FatalIOError);
        }

        for (const label n : nBins_)
        {
            nBin_ *= n;
        }

        if (nBin_ <= 0)
        {
            FatalIOErrorInFunction(binDict)
                << "Number of bins must be greater than zero" << nl
                << "    e1 bins = " << nBins_[0] << nl
                << "    e2 bins = " << nBins_[1] << nl
                << "    e3 bins = " << nBins_[2]
                << exit(FatalIOError);
        }

        Info<< "    - Employing:" << nl
            << "        " << nBins_[0] << " e1 bins," << nl
            << "        " << nBins_[1] << " e2 bins," << nl
            << "        " << nBins_[2] << " e3 bins"
            << endl;

        cumulative_ = binDict.getOrDefault<bool>("cumulative", false);
        Info<< "    - cumulative    : " << cumulative_ << endl;

        if (binDict.found("minMax"))
        {
            const dictionary& minMaxDict(binDict.subDict("minMax"));

            for (direction i = 0; i < vector::nComponents; ++i)
            {
                const word ei("e" + Foam::name(i));

                if (minMaxDict.readIfPresent(ei, binMinMax_[i]))
                {
                    Info<< "    - " << ei << " min        : "
                        << binMinMax_[i][0] << nl
                        << "    - " << ei << " max        : "
                        << binMinMax_[i][1] << endl;
                }
            }
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


void Foam::binModels::uniformBin::initialise()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Determine extents of patches in a given coordinate system
    vector geomMin(GREAT, GREAT, GREAT);
    vector geomMax(-GREAT, -GREAT, -GREAT);

    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = pbm[patchi];
        const vectorField ppcs(csysPtr_->localPosition(pp.faceCentres()));

        for (direction i = 0; i < vector::nComponents; ++i)
        {
            geomMin[i] = min(min(ppcs.component(i)), geomMin[i]);
            geomMax[i] = max(max(ppcs.component(i)), geomMax[i]);
        }
    }

    if (porosity_)
    {
        const HashTable<const porosityModel*> models =
            mesh_.lookupClass<porosityModel>();

        const vectorField dd(csysPtr_->localPosition(mesh_.cellCentres()));

        forAllConstIters(models, iter)
        {
            const porosityModel& pm = *iter();
            const labelList& cellZoneIDs = pm.cellZoneIDs();

            for (const label zonei : cellZoneIDs)
            {
                const cellZone& cZone = mesh_.cellZones()[zonei];
                const vectorField d(dd, cZone);

                for (direction i = 0; i < vector::nComponents; ++i)
                {
                    geomMin[i] = min(min(d.component(i)), geomMin[i]);
                    geomMax[i] = max(max(d.component(i)), geomMax[i]);
                }
            }
        }
    }

    reduce(geomMin, minOp<vector>());
    reduce(geomMax, maxOp<vector>());

    for (direction i = 0; i < vector::nComponents; ++i)
    {
        // Slightly boost max so that region of interest is fully within bounds
        geomMax[i] = 1.0001*(geomMax[i] - geomMin[i]) + geomMin[i];

        // Use geometry limits if not specified by the user
        if (binMinMax_[i][0] == GREAT) binMinMax_[i][0] = geomMin[i];
        if (binMinMax_[i][1] == GREAT) binMinMax_[i][1] = geomMax[i];

        if (binMinMax_[i][0] > binMinMax_[i][1])
        {
            FatalErrorInFunction
                << "Max bounds must be greater than min bounds" << nl
                << "    direction   = " << i << nl
                << "    min         = " << binMinMax_[i][0] << nl
                << "    max         = " << binMinMax_[i][1] << nl
                << exit(FatalError);
        }

        //- Compute bin widths in binning directions
        binW_[i] = (binMinMax_[i][1] - binMinMax_[i][0])/scalar(nBins_[i]);

        if (binW_[i] <= 0)
        {
            FatalErrorInFunction
                << "Bin widths must be greater than zero" << nl
                << "    direction = " << i << nl
                << "    min bound = " << binMinMax_[i][0] << nl
                << "    max bound = " << binMinMax_[i][1] << nl
                << "    bin width = " << binW_[i]
                << exit(FatalError);
        }
    }
}


void Foam::binModels::uniformBin::applyBins
(
    List<Field<vector>>& force,
    List<Field<vector>>& moment,
    const vectorField& d,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT
)
{
    const vectorField dd(csysPtr_->localPosition(d));

    forAll(dd, i)
    {
        // Avoid faces outside of the bin
        bool faceInside = true;
        for (direction j = 0; j < vector::nComponents; ++j)
        {
            if (dd[i][j] < binMinMax_[j][0] || dd[i][j] > binMinMax_[j][1])
            {
                faceInside = false;
                break;
            }
        }

        if (faceInside)
        {
            // Find the bin division corresponding to the face
            Vector<label> n(Zero);
            for (direction j = 0; j < vector::nComponents; ++j)
            {
                n[j] = floor((dd[i][j] - binMinMax_[j][0])/binW_[j]);
                n[j] = min(max(n[j], 0), nBins_[j] - 1);
            }

            // Order: (e1, e2, e3), the first varies the fastest
            const label bini =
                n[0] + nBins_[0]*n[1] + nBins_[0]*nBins_[1]*n[2];

            force[0][bini] += fN[i];
            force[1][bini] += fT[i];
            moment[0][bini] += Md[i]^fN[i];
            moment[1][bini] += Md[i]^fT[i];
        }
    }
}


void Foam::binModels::uniformBin::applyBins
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
    const vectorField dd(csysPtr_->localPosition(d));

    forAll(dd, i)
    {
        // Avoid faces outside of the bin
        bool faceInside = true;
        for (direction j = 0; j < vector::nComponents; ++j)
        {
            if (dd[i][j] < binMinMax_[j][0] || dd[i][j] > binMinMax_[j][1])
            {
                faceInside = false;
                break;
            }
        }

        if (faceInside)
        {
            // Find the bin division corresponding to the face
            Vector<label> n(Zero);
            for (direction j = 0; j < vector::nComponents; ++j)
            {
                n[j] = floor((dd[i][j] - binMinMax_[j][0])/binW_[j]);
                n[j] = min(max(n[j], 0), nBins_[j] - 1);
            }

            // Order: (e1, e2, e3), the first varies the fastest
            const label bini =
                n[0] + nBins_[0]*n[1] + nBins_[0]*nBins_[1]*n[2];

            force[0][bini] += fN[i];
            force[1][bini] += fT[i];
            force[2][bini] += fP[i];
            moment[0][bini] += Md[i]^fN[i];
            moment[1][bini] += Md[i]^fT[i];
            moment[2][bini] += Md[i]^fP[i];
        }
    }
}


void Foam::binModels::uniformBin::createBinnedDataFiles()
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


void Foam::binModels::uniformBin::createBinnedDataFile
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


void Foam::binModels::uniformBin::writeBinnedDataFileHeader
(
    const word& header,
    OFstream& os,
    const bool coeffStyle
) const
{
    writeHeader(os, header + " bins");

    const tensor& R = csysPtr_->R();
    for (direction i = 0; i < vector::nComponents; ++i)
    {
        writeHeaderValue(os, "e"+Foam::name(i)+" bins", nBins_[i]);
        writeHeaderValue(os, "    start", binMinMax_[i][0]);
        writeHeaderValue(os, "    end", binMinMax_[i][1]);
        writeHeaderValue(os, "    delta", binW_[i]);
        writeHeaderValue(os, "    direction", R.col(i));
    }
    writeCommented(os, "bin end co-ords:");
    os  << nl;

    // Compute and print bin end points in binning directions
    for (direction i = 0; i < vector::nComponents; ++i)
    {
        scalar binEnd = binMinMax_[i][0];

        writeCommented(os, "e"+Foam::name(i)+" co-ords   :");
        for (label j = 0; j < nBins_[i]; ++j)
        {
            binEnd += binW_[i];
            os  << tab << binEnd;
        }
        os  << nl;
    }

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


void Foam::binModels::uniformBin::writeBinnedDataFiles
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


void Foam::binModels::uniformBin::writeBinnedDataFile
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


void Foam::binModels::uniformBin::writeBinnedDataFile
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
