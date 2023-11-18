/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

Application
    Test-globalIndex2

Description
    More functional tests for the globalIndex class.

\*---------------------------------------------------------------------------*/

#include "globalIndex.H"
#include "argList.H"
#include "Time.H"
#include "IOstreams.H"
#include "Random.H"
#include "IndirectList.H"
#include "SliceList.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Random rnd(123456);

    // #include "setRootCase.H"
    // #include "createTime.H"

    Info<< nl
        << "simple globalIndex tests" << nl << nl;

    Info<< "Construct from indirect list(s)" << nl;
    {
        // Some sizes
        labelList rawSizes(25);
        forAll(rawSizes, i)
        {
            rawSizes[i] = rnd.position<label>(0, 100);
        }

        Info<< nl
            << "plain sizes: "
            << flatOutput(rawSizes) << nl
            << "    offsets: "
            << flatOutput(globalIndex::calcOffsets(rawSizes))
            << nl;


        sliceRange slice(0, 5, 5);
        Info<< nl
            << "range min/max " << slice.min() << '/' << slice.max() << nl;

        SliceList<label> sliceSizes(rawSizes, slice);

        Info<< nl
            << "indirect addr: " << sliceSizes.addressing() << nl
            << "indirect sizes: "
            << flatOutput(sliceSizes) << nl
            << "       offsets: "
            << flatOutput(globalIndex::calcOffsets(sliceSizes))
            << nl;
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
