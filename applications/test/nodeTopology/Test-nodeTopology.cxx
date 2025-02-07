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

    List<int> interNodeProcs_fake;

    if (UPstream::parRun())
    {
        if (args.found("numProcs"))
        {
            InfoErr<< "ignoring -np option in parallel" << nl;
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

        if (nCores > 1 && nCores < nProcs)
        {
            const label numNodes
                = (nProcs/nCores) + ((nProcs % nCores) ? 1 : 0);

            interNodeProcs_fake.resize(numNodes);

            for (label nodei = 0; nodei < numNodes; ++nodei)
            {
                interNodeProcs_fake[nodei] = nodei * nCores;
            }
        }
    }

    const List<int>& interNodeProcs =
    (
        UPstream::parRun()
      ? UPstream::procID(UPstream::commInterNode())
      : interNodeProcs_fake
    );


    // Generate the graph
    if (UPstream::master(UPstream::worldComm))
    {
        auto& os = Info.stream();

        os << "// node topology graph:" << nl;
        os.beginBlock("graph");

        // Prefer left-to-right layout for large graphs
        os << indent << "rankdir=LR" << nl;

        int pos = 0;

        // First level are the inter-node connections
        const label parent = 0;
        for (const auto proci : interNodeProcs)
        {
            if (parent == proci) continue;

            if (pos)
            {
                os << "  ";
            }
            else
            {
                os << indent;
            }
            os << parent << " -- " << proci;

            if (++pos >= 4)  // Max 4 items per line
            {
                pos = 0;
                os << nl;
            }
        }

        if (pos)
        {
            pos = 0;
            os << nl;
        }

        // Next level are within the nodes
        for (label nodei = 0; nodei < interNodeProcs.size(); ++nodei)
        {
            pos = 0;

            label firstProc = interNodeProcs[nodei];
            const label lastProc =
            (
                (nodei+1 < interNodeProcs.size())
              ? interNodeProcs[nodei+1]
              : nProcs
            );

            os << indent << "// inter-node " << nodei
                << " [" << firstProc
                << ".." << lastProc-1 << "]" << nl;

            for (label proci = firstProc; proci < lastProc; ++proci)
            {
                if (firstProc == proci) continue;

                if (pos)
                {
                    os << "  ";
                }
                else
                {
                    os << indent;
                }
                os << firstProc << " -- " << proci;

                if (++pos >= 4)  // Max 4 items per line
                {
                    pos = 0;
                    os << nl;
                }
            }
            if (pos)
            {
                pos = 0;
                os << nl;
            }
        }

        os.endBlock();
        os << "// end graph" << nl;
    }

    InfoErr << "\nDone" << nl;
    return 0;
}


// ************************************************************************* //
