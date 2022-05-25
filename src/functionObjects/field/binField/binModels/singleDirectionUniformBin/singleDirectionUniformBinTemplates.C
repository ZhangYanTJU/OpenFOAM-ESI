/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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


template<class Type>
void Foam::binModels::singleDirectionUniformBin::writeFileHeader
(
    OFstream& os,
    const bool coeffStyle
) const
{
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
            writeTabbed(os, in + "internal");

            if (decomposePatchValues_)
            {
                writeTabbed(os, in + "normal");
                writeTabbed(os, in + "tangenial");
            }
            else
            {
                writeTabbed(os, in + "patch");
            }
        }
    }
    else
    {
        for (label i = 0; i < nBin_; ++i)
        {
            const word in(Foam::name(i) + ':');
            writeTabbed(os, in + writeComponents<Type>("total"));
            writeTabbed(os, in + writeComponents<Type>("internal"));

            if (decomposePatchValues_)
            {
                writeTabbed(os, in + writeComponents<Type>("normal"));
                writeTabbed(os, in + writeComponents<Type>("tangenial"));
            }
            else
            {
                writeTabbed(os, in + writeComponents<Type>("patch"));
            }
        }
    }

    os  << endl;
}


template<class Type>
bool Foam::binModels::singleDirectionUniformBin::processField
(
    const label fieldi
)
{
    const word& fieldName = fieldNames_[fieldi];

    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const VolFieldType* fieldPtr = mesh_.findObject<VolFieldType>(fieldName);

    if (!fieldPtr)
    {
        return false;
    }

    if (Pstream::master() && !writtenHeader_)
    {
        writeFileHeader<Type>(filePtrs_[fieldi]);
    }

    const VolFieldType& fld = *fieldPtr;

    // Total number of fields
    //
    // 0: internal
    // 1: patch total
    //
    // OR
    //
    // 0: internal
    // 1: patch normal
    // 2: patch tangential
    label nField = 2;
    if (decomposePatchValues_)
    {
        nField += 1;
    }

    List<List<Type>> data(nField);
    for (auto& binList : data)
    {
        binList.setSize(nBin_, Zero);
    }

    for (const label zonei : cellZoneIDs_)
    {
        const cellZone& cZone = mesh_.cellZones()[zonei];

        for (const label celli : cZone)
        {
            scalar dd = mesh_.C()[celli] & binDir_;

            if (dd < binMin_ || dd > binMax_)
            {
                continue;
            }

            // Find the bin division corresponding to the cell
            const label bini =
                min(max(floor((dd - binMin_)/binDx_), 0), nBin_ - 1);

            data[0][bini] += fld[celli];
        }
    }

    forAllIters(patchSet_, iter)
    {
        const label patchi = iter();
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const vectorField np(mesh_.boundary()[patchi].nf());

        const scalarField dd(pp.faceCentres() & binDir_);

        forAll(dd, facei)
        {
            // Avoid faces outside of the bin
            if (dd[facei] < binMin_ || dd[facei] > binMax_)
            {
                continue;
            }

            // Find the bin division corresponding to the face
            const label bini =
                min(max(floor((dd[facei] - binMin_)/binDx_), 0), nBin_ - 1);

            const Type& v = fld.boundaryField()[patchi][facei];

            if (!decomposePatchValues(data, bini, v, np[facei]))
            {
                data[1][bini] += v;
            }
        }
    }

    if (Pstream::master())
    {
        writeBinnedData(data, filePtrs_[fieldi]);
    }

    return true;
}


template<class Type>
void Foam::binModels::singleDirectionUniformBin::writeBinnedData
(
    List<List<Type>>& data,
    Ostream& os
) const
{
    if (cumulative_)
    {
        for (auto& datai : data)
        {
            for (label bini = 1; bini < nBin_; ++bini)
            {
                datai[bini] += datai[bini-1];
            }
        }
    }

    writeCurrentTime(os);

    for (label bini = 0; bini < nBin_; ++bini)
    {
        Type total = Zero;

        for (label i = 0; i < data.size(); ++i)
        {
            total += data[i][bini];
        }

        os  << tab << total;

        for (label i = 0; i < data.size(); ++i)
        {
            os  << tab << data[i][bini];
        }
    }

    os  << nl;
}


// ************************************************************************* //
