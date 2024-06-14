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
    createViewFactors

Group
    grpPreProcessingUtilities

Description
    Creates view factors for the view factor radiation model.

    User-selectable models:

    - raySearchEngine: model to generate rays, i.e. face-to-face connections
    - viewFactorModel: model to compute the view factors

    For visualisation, use:

    - Write the view factors as a volume field

        writeViewFactors    yes;

    - Write the rays using OBJ format:

        writeRays       yes; // default = no

    Participating patches must be in the \c vewFactorWall group, i.e. using the
    \c inGroups entry of the "\<case\>/polyMesh/boundary" file.

    \verbatim
    myPatch
    {
        type            wall;
        inGroups        2(wall viewFactorWall);
        ...
    }
    \endverbatim

    Reads:

    - <constant>/viewFactorsDict : main controls
    - <constant>/finalAgglom : agglomeration addressing (from faceAgglomerate)


    Generates:

    - <constant>/F : view factors (matrix)
    - <constant>/mapDist : map used for parallel running
    - <constant>/globalFaceFaces : face addressing

SeeAlso
- Foam::VF::raySearchEngine
- Foam::VE::viewFactorModel

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "viewFactorModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "viewFactorsDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    );


    // Calculate the view factors
    auto modelPtr = VF::viewFactorModel::New(mesh, dict);

    modelPtr->calculate();

    Info<< nl;

    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //