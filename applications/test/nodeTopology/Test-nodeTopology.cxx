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

Application
    Test-nodeTopology

Description
    Simple reporting of node topology

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();
    argList::addOption
    (
        "numProcs",
        "int",
        "Num of ranks to simulate (default: 16)"
    );
    argList::addOption
    (
        "cores",
        "int",
        "Num of cores to simulate (default: 4)"
    );

    #include "setRootCase.H"

    label nProcs = UPstream::nProcs(UPstream::worldComm);

    DynamicList<int> fake_interNode_offsets;

    if (UPstream::parRun())
    {
        if (args.found("numProcs"))
        {
            InfoErr<< "ignoring -numProcs option in parallel" << nl;
        }
        if (args.found("cores"))
        {
            InfoErr<< "ignoring -cores option in parallel" << nl;
        }
    }
    else
    {
        // serial
        nProcs = args.getOrDefault<label>("numProcs", 16);
        label nCores = args.getOrDefault<label>("cores", 4);

        auto& interNode_offsets = fake_interNode_offsets;

        if (nCores > 1 && nCores < nProcs)
        {
            // Build the inter-node offsets
            interNode_offsets.reserve((nProcs/nCores) + 4);
            interNode_offsets.push_back(0);

            for
            (
                int count = interNode_offsets.back() + nCores;
                count < nProcs;
                count += nCores
            )
            {
                interNode_offsets.push_back(count);
            }

            interNode_offsets.push_back(nProcs);
        }
        else
        {
            // Some fallback
            interNode_offsets.reserve(2);
            interNode_offsets.push_back(0);
            interNode_offsets.push_back(nProcs);
        }
    }

    const List<int>& interNodeOffsets =
    (
        UPstream::parRun()
      ? UPstream::interNode_offsets()
      : fake_interNode_offsets
    );


    if (UPstream::parRun())
    {
        const auto& procs = UPstream::localNode_parentProcs();
        Perr<< "local processors: [" << procs.min()
            << ".." << procs.max() << ']' << endl;
    }

    // Generate the graph
    if (UPstream::master(UPstream::worldComm))
    {
        auto& os = Info.stream();

        os << "// node topology graph:" << nl;
        os.beginBlock("graph");

        // Prefer left-to-right layout for large graphs
        os << indent << "rankdir=LR" << nl;

        const label numNodes = interNodeOffsets.size()-1;

        // First level are the inter-node connections
        {
            os  << indent << 0 << " -- " << token::LBRACE;

            for (label nodei = 1; nodei < numNodes; ++nodei)
            {
                os  << ' ' << interNodeOffsets[nodei];
            }

            os  << token::SPACE << token::RBRACE
                << "  // inter-node: " << flatOutput(interNodeOffsets)
                << nl;
        }

        // Next level are the local-node connections
        for (label nodei = 0; nodei < numNodes; ++nodei)
        {
            const auto firstProc = interNodeOffsets[nodei];
            const auto lastProc = interNodeOffsets[nodei+1];

            os  << indent << firstProc << " -- " << token::DQUOTE
                << (firstProc+1) << ".." << (lastProc-1)
                << token::DQUOTE << nl;
        }

        os.endBlock();
        os << "// end graph" << nl;
    }

    InfoErr << "\nDone" << nl;
    return 0;
}


// ************************************************************************* //
