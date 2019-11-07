/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Test-SubField

Description
    Simple tests on SubList, SubField

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "argList.H"

#include "scalarField.H"
#include "SubField.H"
#include "labelRange.H"
#include <numeric>

using namespace Foam;

template<class T>
void print(const UList<T>& list)
{
    Info<< flatOutput(list) << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();

    {
        List<scalar> ident(25);
        std::iota(ident.begin(), ident.end(), 0);

        print(ident);

        SubList<scalar>(ident, 10) = -10;
        print(ident);

        SubField<scalar>(ident, 10) = 10;
        print(ident);

        SubField<scalar>(ident, 10) += 10;
        print(ident);

        SubField<scalar>{ident, 10, 10} *= 5;
        print(ident);


        // NOTE: Need {} instead of ()
        // SubList<scalar>(ident) = 100;

        // GCC
        // error: conflicting declaration 'Foam::SubList<double> ident'

        // CLANG
        // warning: parentheses were disambiguated as redundant parentheses
        // around declaration of variable named 'ident' [-Wvexing-parse]

        SubList<scalar>{ident} = 100;
        print(ident);
    }

    Info << "\nEnd\n";
    return 0;
}

// ************************************************************************* //