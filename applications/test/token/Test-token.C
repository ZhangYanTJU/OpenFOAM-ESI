/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Description
    Test token construct assign etc.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "StringStream.H"
#include "cpuTime.H"
#include "labelList.H"
#include "scalarList.H"
#include "DynamicList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();

    argList args(argc, argv, false, true);

    token tok1;
    Info<< "default construct: " << tok1.info() << endl;

    tok1 = double(3.14159);
    Info<< "assign double: " << tok1.info() << endl;

    token tok2(tok1);
    Info<< "copy construct: " << tok2.info() << endl;

    tok1 = word("this-word");
    Info<< "assign word: " << tok1.info() << endl;

    token tok3(tok1);
    Info<< "copy construct: " << tok3.info() << endl;
    Info<< "orig: " << tok1.info() << endl;

    token tok4(std::move(tok1));
    Info<< "move construct: " << tok4.info() << endl;
    Info<< "orig: " << tok1.info() << endl;

    tok3 = tok4;
    Info<< "assign token: " << tok3.info() << endl;
    Info<< "orig: " << tok4.info() << endl;

    //
    // Compound
    //

    {
        // This version is good

        token ctok1(new token::Compound<labelList>(identity(10)));

        Info<< "compound from pointer: "
            << ctok1.info() << nl << ctok1 << endl;
    }

    {
        // This also works, but not actually using the autoPtr directly

        autoPtr<token::Compound<labelList>> ptr
        (
            new token::Compound<labelList>(identity(10, -9))
        );

        token ctok1(ptr.release());  // release() not get()!

        Info<< "compound from autoPtr: "
            << ctok1.info() << nl << ctok1 << endl;
    }

    {
        // Construct from pointer
        autoPtr<token::compound> ptr
        (
            token::compound::New("List<label>")
        );

        token ctok1(ptr.release());  // release() not get()!

        Info<< "compound from New (via pointer): "
            << ctok1.info() << nl << ctok1 << endl;
    }

    {
        // Construct from autoPtr
        autoPtr<token::compound> ptr
        (
            token::Compound<scalarList>::New(10, 1.0)
        );

        token ctok1(std::move(ptr));
        Info<< "compound from autoPtr: "
            << ctok1.info() << nl << ctok1 << endl;

        // Shrink
        ctok1.refCompoundToken().resize(5);

        Info<< "resized: "
            << ctok1.info() << nl << ctok1 << endl;

        const scalarList* listptr = ctok1.compoundToken().isA<scalarList>();
        if (listptr)
        {
            for (scalar& val : const_cast<scalarList&>(*listptr))
            {
                val *= 5;
            }

            Info<< "multiplied List<scalar>: "
                << ctok1.info() << nl << ctok1 << endl;
        }

        listptr = ctok1.isCompound<scalarList>();
        if (listptr)
        {
            for (scalar& val : const_cast<scalarList&>(*listptr))
            {
                val /= 2;
            }

            Info<< "divided List<scalar>: "
                << ctok1.info() << nl << ctok1 << endl;
        }

        const labelList* listptr2 = ctok1.isCompound<labelList>();
        if (listptr2)
        {
            for (label& val : const_cast<labelList&>(*listptr2))
            {
                val /= 2;
            }

            Info<< "divided List<label>: "
                << ctok1.info() << nl << ctok1 << endl;
        }
        else
        {
            Info<< "compound is not List<label>" << nl;
        }

        Info<< "Before fill_zero: " << ctok1 << endl;

        ctok1.refCompoundToken().fill_zero();

        Info<< "After fill_zero: " << ctok1 << endl;


        if (ctok1.isCompound())
        {
            auto& ct = ctok1.refCompoundToken();

            ct.resize(20);
            bool handled = true;

            switch (ct.typeCode())
            {
                case token::tokenType::BOOL :
                {
                    UList<bool> cmpts
                    (
                        reinterpret_cast<bool*>(ct.data_bytes()),
                        label(ct.size_bytes() / sizeof(bool))
                    );
                    cmpts = false;
                }
                break;

                case token::tokenType::LABEL :
                {
                    UList<label> cmpts
                    (
                        reinterpret_cast<label*>(ct.data_bytes()),
                        label(ct.size_bytes() / sizeof(label))
                    );
                    cmpts = 123;
                }
                break;

                case token::tokenType::FLOAT :
                {
                    UList<float> cmpts
                    (
                        reinterpret_cast<float*>(ct.data_bytes()),
                        label(ct.size_bytes() / sizeof(float))
                    );
                    cmpts = 2.7;
                }
                break;

                case token::tokenType::DOUBLE :
                {
                    UList<double> cmpts
                    (
                        reinterpret_cast<double*>(ct.data_bytes()),
                        label(ct.size_bytes() / sizeof(double))
                    );
                    cmpts = 3.1415;
                }
                break;

                default:
                    handled = false;
                    break;
            }


            if (handled)
            {
                Info<< "assigned: " << ctok1 << nl;
            }
        }
    }

    return 0;
}


// ************************************************************************* //
