/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
    Test-CompactListList

Description
    Simple demonstration and test application for the CompactListList class.

\*---------------------------------------------------------------------------*/

#include "CompactListList.H"
#include "IndirectList.H"
#include "IOstreams.H"
#include "SpanStream.H"
#include "faceList.H"

#include <iterator>  // for back_inserter

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    {
        // Default construct
        CompactListList<label> cll1;
        Info<< "cll1:" << cll1 << nl;

        // Resize and assign row by row
        labelList row0(2, Zero);
        labelList row1(3, label(1));

        labelList rowSizes(2);
        rowSizes[0] = row0.size();
        rowSizes[1] = row1.size();
        cll1.resize(rowSizes);

        cll1[0].deepCopy(row0);
        cll1[1].deepCopy(row1);
        Info<< "cll1:" << cll1 << endl;

        forAll(cll1.values(), i)
        {
            Info<< "i:" << i << " whichRow:" << cll1.whichRow(i) << endl;
        }

        Info<< "unpack:" << cll1.unpack<face>() << endl;
    }

    List<List<label>> lll(5);
    lll[0].setSize(3, 0);
    lll[1].setSize(2, 1);
    lll[2].setSize(6, 2);
    lll[3].setSize(0, 3);
    lll[4].setSize(1, 4);

    Info<< "packed:" << CompactListList<label>::pack(lll) << endl;

    auto cll2(CompactListList<label>::pack(lll));

    Info<< "cll2  = " << cll2 << endl;

    forAll(cll2, i)
    {
        Info<< cll2[i] << nl;
    }
    Info<< endl;

    Info<< "cll2(2, 3) = " << cll2(2, 3) << nl << endl;
    cll2(2, 3) = 999;
    Info<< "cll2(2, 3) = " << cll2(2, 3) << nl << endl;

    Info<< "cll2 as List<List<label>> " << cll2.unpack() << endl;

    cll2.setSize(3);

    Info<< "cll2  = " << cll2 << endl;

    cll2.setSize(0);

    Info<< "cll2  = " << cll2 << endl;


    List<label> rowSizes(5);
    rowSizes[0] = 2;
    rowSizes[1] = 0;
    rowSizes[2] = 1;
    rowSizes[3] = 3;
    rowSizes[4] = 2;

    CompactListList<label> cll3(rowSizes, 1);

    Info<< "cll3 = " << cll3 << endl;

    CompactListList<label> cll4;

    cll4.transfer(cll3);

    Info<< "cll3 = " << cll3 << endl;
    Info<< "cll4 = " << cll4 << endl;


    {
        // IO
        OCharStream ostr;
        ostr << cll4;

        ISpanStream istr(ostr.view());
        CompactListList<label> cll5(istr);
        Info<< "cll5 = " << cll5 << endl;
    }
    {
        // IO
        cll4.clear();
        OCharStream ostr;
        ostr << cll4;

        ISpanStream istr(ostr.view());
        CompactListList<label> cll5(istr);
        Info<< "cll5 = " << cll5 << endl;
    }

    // Make some faces
    {
        faceList fcs(5);
        forAll(fcs, facei)
        {
            fcs[facei] = face(identity(4, facei));
        }

        Info<< "input faces: " << fcs << endl;

        // From <face>
        auto compactFcs(CompactListList<label>::pack<face>(fcs));
        Info<< "compact faces:" << compactFcs << endl;

        faceList fcs2 = compactFcs.unpack<face>();
        Info<< "deserialized:" << fcs2 << endl;

        // Unpack some faces
        DynamicList<face> extracted(compactFcs.size());

        compactFcs.copy_unpack<face>
        (
            std::back_inserter(extracted),
            2, 2
        );

        Info<< "copy_unpack 1: " << extracted << nl;

        compactFcs.copy_unpack<face>
        (
            std::back_inserter(extracted)
            // labelRange(2, 1)
        );

        Info<< "copy_unpack 2: " << extracted << nl;

        // From some faces
        IndirectList<face> subfaces(fcs, labelList({2, 4, 1}));

        Info<< "sub faces: " << subfaces << endl;

        auto subCompact
        (
            CompactListList<label>::pack(subfaces)
        );
        Info<< "compact faces:" << subCompact << endl;
        Info<< "matrix content:" << nl;
        subCompact.writeMatrix(Info) << endl;
    }

    return 0;
}


// ************************************************************************* //
