/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-exprValue

Description
    Test low-level polymorphic value container (exprValue)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "ITstream.H"
#include "exprValue.H"

using namespace Foam;

void printInfo(const expressions::exprValue& val)
{
    Info<< "Boxed type:" << int(val.typeCode())
        << " (" << val.valueTypeName() << ") good:"
        << val.good() << " => " << val << nl;
}


expressions::exprValue tryParse(const std::string& str)
{
    expressions::exprValue val, val2;

    ITstream is(str);

    const bool ok = val.read(is);

    Info<< "read " << Foam::name(val.typeCode()) << " from " << str;

    if (ok)
    {
        Info<< " trailing tokens:" << is.nRemainingTokens() << nl
            << "value: " << val << nl;
    }
    else
    {
        Info<< " FAILED" << nl;
    }

    if (ok)
    {
        Info<< "Direct from string: ";
        if (expressions::exprValue::read(str, val2))
        {
            Info<< "good" << nl;
        }
        else
        {
            Info<< "bad" << nl;
        }
    }
    return val;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();

    #include "setRootCase.H"

    // Aborts
    // expressions::exprValue value(std::string(""));

    {
        expressions::exprValue value;

        // Nothing
        printInfo(value);

        value.set(scalar(100));
        printInfo(value);

        value.set(vector(1,2,3));
        printInfo(value);

        value = vector(4,5,6);
        printInfo(value);

        value = Zero;
        printInfo(value);

        value.clear();
        printInfo(value);

        value = 100 * vector(1,0,0);
        printInfo(value);
    }


    {
        Info<< nl << "Test parsing" << nl << nl;

        for
        (
            const auto& input :
            stringList
            ({
                "()",  // bad
                "( 1 2  ",  // also bad
                "(  ",  // really bad
                "(1 16 12)",
                "(1 bad)",
                "(5)",
                "1.2345",
                "5.678 trailing",
                "true",
                "false",
                " 1 ",
                "  yes no "
            })
        )
        {
            (void) tryParse(input);
        }
    }

    return 0;
}


// ************************************************************************* //
