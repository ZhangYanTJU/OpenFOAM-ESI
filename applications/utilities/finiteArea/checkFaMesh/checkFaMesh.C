/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd.
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

Application
    checkFaMesh

Description
    Inspects and reports mesh metrics for a given finiteArea mesh.

Usage
    \b checkFaMesh [OPTIONS]

    \param -region \<name\> \n
        Specify an alternative mesh region.

    The inherited options are elaborated in argList.H.

Note
    - \c checkFaMesh does not check the validity
    of a given finiteArea mesh unlike \c checkMesh.
    - Parallel execution is currently not available.

See also
    Foam::makeFaMesh

Author
    Zeljko Tukovic, FAMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Inspects and reports mesh metrics for a given finiteArea mesh."
    );

    #include "addRegionOption.H"
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"
    #include "createFaMesh.H"

    Info<< "Time = " << runTime.timeName() << nl << endl;

    Info<< "Mesh stats:" << nl
        << "    points: " << aMesh.nPoints() << nl
        << "    internal edges: " << aMesh.nInternalEdges() << nl
        << "    edges: " << aMesh.nEdges() << nl
        << "    faces: " << aMesh.nFaces() << nl
        << endl;

    Info<< "Geometry stats:" << nl
        << "    Face area:"
        << " min = " << min(aMesh.S().field())
        << " max = " << max(aMesh.S().field()) << nl
        << "    Internal edge length:"
        << " min = " << min(aMesh.magLe().internalField()).value()
        << " max = " << max(aMesh.magLe().internalField()).value() << nl
        << "    Edge length:"
        << " min = " << min(aMesh.magLe()).value()
        << " max = " << max(aMesh.magLe()).value() << nl
        << "    Face area normals:"
        << " min = " << min(aMesh.faceAreaNormals().field())
        << " max = " << max(aMesh.faceAreaNormals().field()) << nl
        << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
