/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "noBin.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{
    defineTypeNameAndDebug(noBin, 0);
    addToRunTimeSelectionTable(binModel, noBin, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binModels::noBin::noBin
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    binModel(name, dict, mesh)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::binModels::noBin::read(const dictionary& dict)
{
    return binModel::read(dict);
}


void Foam::binModels::noBin::initialise()
{
    nBin_ = 1;
}


void Foam::binModels::noBin::applyBins
(
    List<Field<vector>>& force,
    List<Field<vector>>& moment,
    const vectorField& d,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT
)
{
    force[0][0] += sum(fN);
    force[1][0] += sum(fT);
    moment[0][0] += sum(Md^fN);
    moment[1][0] += sum(Md^fT);
}


void Foam::binModels::noBin::applyBins
(
    List<Field<vector>>& force,
    List<Field<vector>>& moment,
    const vectorField& d,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT,
    const vectorField& fP
)
{
    force[0][0] += sum(fN);
    force[1][0] += sum(fT);
    force[2][0] += sum(fP);
    moment[0][0] += sum(Md^fN);
    moment[1][0] += sum(Md^fT);
    moment[2][0] += sum(Md^fP);
}


// ************************************************************************* //
