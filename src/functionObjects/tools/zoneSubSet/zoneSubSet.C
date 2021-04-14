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

#include "zoneSubSet.H"
#include "topoSetSource.H"
#include "topoSet.H"
#include "haloToCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoneSubSet, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneSubSet::init(const dictionary& dict)
{
    dict.readIfPresent("cellZones", zoneNames_);

    if (zoneNames_.size())
    {
        autoPtr<topoSet> currentSet;
        const topoSetSource::setAction action = topoSetSource::ADD;
        const word setType("cellZoneSet");

        currentSet = topoSet::New
        (
            setType,
            mesh_,
            zoneNames_[0],
            IOobject::MUST_READ
        );

        if (debug)
        {
            Pout<< "Read set " << currentSet().type() << ' '
                << zoneNames_[0] << " with size "
                << currentSet().size()
                << endl;
        }

        dictionary dictHalo;
        dictHalo.add("steps", nLayers_);

        autoPtr<topoSetSource> haloSource = topoSetSource::New
        (
            "haloToCell",
            mesh_,
            dictHalo
        );

        haloSource().verbose(false);
        haloSource().applyToSet(action, currentSet());

        for (label zonei = 1; zonei < zoneNames_.size(); ++zonei)
        {
            dictionary dict;
            dict.add("zone", zoneNames_[zonei]);

            autoPtr<topoSetSource> source = topoSetSource::New
            (
                "zoneToCell",
                mesh_,
                dict
            );
            source().applyToSet(action, currentSet());

            dictionary dictHalo;
            dictHalo.add("steps", nLayers_);
            haloSource.reset
            (
                topoSetSource::New
                (
                    "haloToCell",
                    mesh_,
                    dictHalo
                )
            );
            haloSource().verbose(false);
            haloSource().applyToSet(action, currentSet());
        }

        Info<< "Created " << currentSet().type() << ' ' << setType << endl;

        if (currentSet)
        {
            Info<< "    "
                << currentSet().type() << ' '
                << currentSet().name() << " now size "
                << returnReduce(currentSet().size(), sumOp<label>())
                << endl;
        }

        subSetMeshPtr_.reset
        (
            new fvMeshSubset
            (
                mesh_,
                currentSet()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneSubSet::zoneSubSet
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    nLayers_(dict.getCheckOrDefault<label>("nLayers", 0, labelMinMax::ge(0)))
{
    init(dict);
}


// ************************************************************************* //
