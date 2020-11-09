/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "argListRedirect.H"
#include "IOstreams.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

inline bool opt_join(const char* optName)
{
    return strcmp(optName, "join-stderr") == 0;
}

inline bool opt_rank(const char* optName)
{
    return strcmp(optName, "append-rank") == 0;
}

inline bool opt_stderr(const char* optName)
{
    return strcmp(optName, "stderr") == 0;
}

inline bool opt_stdout(const char* optName)
{
    return strcmp(optName, "stdout") == 0;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Detail::redirectOutputs::redirectOutputs(int& argc, char**& argv)
:
    stdout_(),
    stderr_(),
    join_(false),
    ranks_(false)
{
    List<bool> skip(label(argc), false);
    bool filter = false;

    for (int argi = 1; argi < argc-1; ++argi)
    {
        if (argv[argi][0] == '-')
        {
            const char *optName = &argv[argi][1];

            if (opt_join(optName))
            {
                join_ = true;
                filter = true;
                skip[argi] = true;
            }
            else if (opt_rank(optName))
            {
                ranks_ = true;
                filter = true;
                skip[argi] = true;
            }
            else if (opt_stdout(optName))
            {
                stdout_ = argv[argi+1];
                filter = true;
                skip[argi] = true;
                skip[argi+1] = true;
                ++argi;
            }
            else if (opt_stderr(optName))
            {
                stderr_ = argv[argi+1];
                filter = true;
                skip[argi] = true;
                skip[argi+1] = true;
                ++argi;
            }
        }
    }


    // Test final arg separately
    {
        const int argi = argc-1;

        if (argi > 0 && argv[argi][0] == '-')
        {
            const char *optName = &argv[argi][1];

            if (opt_join(optName))
            {
                join_ = true;
                filter = true;
                skip[argi] = true;
            }
            else if (opt_rank(optName))
            {
                ranks_ = true;
                filter = true;
                skip[argi] = true;
            }
        }
    }


    if (filter)
    {
        int nArgs = 1;

        for (int argi = 1; argi < argc; ++argi)
        {
            if (!skip[argi])
            {
                argv[nArgs] = argv[argi];
                ++nArgs;
            }
        }
        argc = nArgs;
    }


    // Resolve potential conflicts

    if (!stderr_.empty())
    {
        if (stdout_.empty())
        {
            join_ = false;
        }
        else if (stdout_ == stderr_)
        {
            join_ = true;
            stderr_.clear();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::Detail::redirectOutputs::active() const
{
    return join_ || !stdout_.empty() || !stderr_.empty();
}


// ************************************************************************* //
