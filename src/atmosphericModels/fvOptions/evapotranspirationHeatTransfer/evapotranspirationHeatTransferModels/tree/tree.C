/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "tree.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace evapotranspirationHeatTransferModels
{
    defineTypeNameAndDebug(tree, 0);
    addToRunTimeSelectionTable
    (
        evapotranspirationHeatTransferModel,
        tree,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::evapotranspirationHeatTransferModels::tree::Et() const
{
    return a_*q() + b_;
}


Foam::tmp<Foam::scalarField>
Foam::evapotranspirationHeatTransferModels::tree::leafArea
(
    const labelList& cells
) const
{
    const volScalarField& LAD = getOrReadField(LADname_);

    const scalarField& V = mesh().V();

    auto tleafArea = tmp<scalarField>::New(cells.size(), Zero);
    auto& leafArea = tleafArea.ref();

    forAll(cells, i)
    {
        const label celli = cells[i];

        leafArea[i] = LAD[celli]*V[celli];
    }

    return tleafArea;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evapotranspirationHeatTransferModels::tree::tree
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    evapotranspirationHeatTransferModel(dict, mesh),
    a_(),
    b_(),
    lambda_(),
    LADname_()
{
    Info<< "    Activating evapotranspiration heat transfer model: "
        << typeName << endl;

    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::evapotranspirationHeatTransferModels::tree::Q
(
    const labelList& cells
) const
{
    // Convert units from [MJ g/hr] to [J kg/s]
    static const scalar unitConverter = scalar(1000)/scalar(3600);

    return unitConverter*lambda_*Et()*leafArea(cells);
}


bool Foam::evapotranspirationHeatTransferModels::tree::read
(
    const dictionary& dict
)
{
    if (!evapotranspirationHeatTransferModel::read(dict))
    {
        return false;
    }

    const auto range = scalarMinMax::ge(SMALL);

    a_ = dict.getOrDefault<scalar>("a", 0.3622);
    b_ = dict.getOrDefault<scalar>("b", 60.758);
    lambda_ = dict.getCheckOrDefault<scalar>("lambda", 2.44, range);
    LADname_ = dict.getOrDefault<word>("LAD", "LAD");

    (void) getOrReadField(LADname_);

    return true;
}


// ************************************************************************* //
