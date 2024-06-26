/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "sizeDistribution.H"
#include "sizeGroup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sizeDistribution, 0);
    addToRunTimeSelectionTable(functionObject, sizeDistribution, dictionary);
}
}


const Foam::Enum
<
    Foam::functionObjects::sizeDistribution::selectionModeTypes
>
Foam::functionObjects::sizeDistribution::selectionModeTypeNames_
({
    {selectionModeTypes::rtCellZone, "cellZone"},
    {selectionModeTypes::rtAll, "all"},
});


const Foam::Enum
<
    Foam::functionObjects::sizeDistribution::functionTypes
>
Foam::functionObjects::sizeDistribution::functionTypeNames_
({
    {functionTypes::ftNdf, "numberDensity"},
    {functionTypes::ftVdf, "volumeDensity"},
    {functionTypes::ftNc, "numberConcentration"},
    {functionTypes::ftMom, "moments"},
});


const Foam::Enum
<
    Foam::functionObjects::sizeDistribution::abszissaTypes
>
Foam::functionObjects::sizeDistribution::abszissaTypeNames_
({

    {abszissaTypes::atDiameter, "diameter"},
    {abszissaTypes::atVolume, "volume"},
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::sizeDistribution::initialise
(
    const dictionary& dict
)
{
    switch (functionType_)
    {
        case ftNdf:
        {
            break;
        }

        case ftVdf:
        {
            break;
        }

        case ftNc:
        {
            break;
        }

        case ftMom:
        {
            break;
        }

        default:
        {
            FatalIOErrorInLookup
            (
                dict,
                "functionType",
                word::null,
                functionTypeNames_
            ) << exit(FatalIOError);
        }
    }

    switch (abszissaType_)
    {
        case atDiameter:
        {
            break;
        }

        case atVolume:
        {
            break;
        }

        default:
        {
            FatalIOErrorInLookup
            (
                dict,
                "abszissaType",
                word::null,
                abszissaTypeNames_
            ) << exit(FatalIOError);
        }
    }

    setCellZoneCells();

    if (nCells_ == 0)
    {
        FatalIOErrorInFunction(dict)
            << type() << " " << name() << ": "
            << selectionModeTypeNames_[selectionModeType_]
            << '(' << selectionModeTypeName_ << "):" << nl
            << "    Selection has no cells" << exit(FatalIOError);
    }

    volume_ = volume();

    Info<< type() << " " << name() << ":"
        << selectionModeTypeNames_[selectionModeType_]
        << '(' << selectionModeTypeName_ << "):" << nl
        << "    total cells  = " << nCells_ << nl
        << "    total volume = " << volume_
        << nl << endl;
}


void Foam::functionObjects::sizeDistribution::setCellZoneCells()
{
    switch (selectionModeType_)
    {
        case rtCellZone:
        {
            dict().readEntry("cellZone", selectionModeTypeName_);

            label zoneId =
                mesh().cellZones().findZoneID(selectionModeTypeName_);

            if (zoneId < 0)
            {
                FatalIOErrorInFunction(dict_)
                    << "Unknown cellZone name: " << selectionModeTypeName_
                    << ". Valid cellZone names are: "
                    << mesh().cellZones().names()
                    << nl << exit(FatalIOError);
            }

            cellId_ = mesh().cellZones()[zoneId];
            nCells_ = returnReduce(cellId_.size(), sumOp<label>());
            break;
        }

        case rtAll:
        {
            cellId_ = identity(mesh().nCells());
            nCells_ = returnReduce(cellId_.size(), sumOp<label>());
            break;
        }

        default:
        {
            FatalIOErrorInLookup
            (
                dict_,
                "selectionMode",
                word::null,
                selectionModeTypeNames_
            ) << exit(FatalIOError);
        }
    }
}


Foam::scalar Foam::functionObjects::sizeDistribution::volume() const
{
    return gSum(filterField(mesh().V()));
}


void Foam::functionObjects::sizeDistribution::combineFields(scalarField& field)
{
    List<scalarField> allValues(Pstream::nProcs());

    allValues[Pstream::myProcNo()] = field;

    Pstream::gatherList(allValues);

    if (Pstream::master())
    {
        field =
            ListListOps::combine<scalarField>
            (
                allValues,
                accessOp<scalarField>()
            );
    }
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::sizeDistribution::filterField
(
    const scalarField& field
) const
{
    return tmp<scalarField>::New(field, cellId_);
}


void Foam::functionObjects::sizeDistribution::writeFileHeader
(
    const label i
)
{
    OFstream& file = this->file();

    switch (functionType_)
    {
        case ftNdf:
        {
            writeHeader(file, "Number density function");
            break;
        }

        case ftVdf:
        {
            writeHeader(file, "Volume density function");
            break;
        }

        case ftNc:
        {
            writeHeader(file, "Number concentration");
            break;
        }

        case ftMom:
        {
            writeHeader(file, "Moments");
            break;
        }
    }

    switch (abszissaType_)
    {
        case atVolume:
        {
            writeCommented(file, "Time/volume");
            break;
        }

        case atDiameter:
        {
            writeCommented(file, "Time/diameter");
            break;
        }
    }

    switch (functionType_)
    {
        case ftMom:
        {
            for (label i = 0; i <= momentOrder_; i++)
            {
                file() << tab << i;
            }

            break;
        }
        default:
        {
            forAll(popBal_.sizeGroups(), sizeGroupi)
            {
                const diameterModels::sizeGroup& fi =
                    popBal_.sizeGroups()[sizeGroupi];

                switch (abszissaType_)
                {
                    case atDiameter:
                    {
                        file() << tab  << fi.d().value();

                        break;
                    }

                    case atVolume:
                    {
                        file() << tab  << fi.x().value();

                        break;
                    }
                }
            }

            break;
        }
    }

    file << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sizeDistribution::sizeDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name),
    dict_(dict),
    selectionModeType_
    (
        selectionModeTypeNames_.get("selectionMode", dict)
    ),
    selectionModeTypeName_(),
    functionType_(functionTypeNames_.get("functionType", dict)),
    abszissaType_(abszissaTypeNames_.get("abszissaType", dict)),
    nCells_(0),
    cellId_(),
    volume_(0.0),
    writeVolume_(dict.getOrDefault("writeVolume", false)),
    popBal_
    (
        obr_.lookupObject<Foam::diameterModels::populationBalanceModel>
        (
            dict.get<word>("populationBalance")
        )
    ),
    N_(popBal_.sizeGroups().size()),
    momentOrder_(dict.getOrDefault<label>("momentOrder", 0)),
    normalize_(dict.getOrDefault("normalize", false)),
    sumN_(0.0),
    sumV_(0.0)
{
    read(dict);
    resetFile(name);
    newFileAtStartTime(name);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sizeDistribution::~sizeDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sizeDistribution::read(const dictionary& dict)
{
    if (dict != dict_)
    {
        dict_ = dict;
    }

    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    initialise(dict);

    return true;
}


bool Foam::functionObjects::sizeDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::sizeDistribution::write()
{
    writeFileHeader();
    writeCurrentTime(file());

    Log << type() << " " << name() << " write" << nl;

    scalarField V(filterField(mesh().V()));
    combineFields(V);

    sumN_ = 0;
    sumV_ = 0;

    forAll(N_, i)
    {
        const Foam::diameterModels::sizeGroup& fi = popBal_.sizeGroups()[i];

        const volScalarField& alpha = fi.VelocityGroup().phase();

        scalarField Ni(fi*alpha/fi.x());
        scalarField values(filterField(Ni));
        scalarField V(filterField(mesh().V()));

        // Combine onto master
        combineFields(values);
        combineFields(V);

        if (Pstream::master())
        {
            // Calculate volume-averaged number concentration
            N_[i] = sum(V*values)/sum(V);
        }

        sumN_ += N_[i];

        sumV_ += N_[i]*fi.x().value();
    }

    if (Pstream::master())
    {
        switch (functionType_)
        {
            case ftMom:
            {
                for (label m = 0; m <= momentOrder_; m++)
                {
                    scalar result(0.0);

                    forAll(N_, i)
                    {
                        const Foam::diameterModels::sizeGroup& fi =
                            popBal_.sizeGroups()[i];

                        switch (abszissaType_)
                        {
                            case atVolume:
                            {
                                result += pow(fi.x().value(), m)*N_[i];

                                break;
                            }

                            case atDiameter:
                            {
                                result += pow(fi.d().value(), m)*N_[i];

                                break;
                            }
                        }
                    }

                    file() << tab << result;
                }

                break;
            }

            default:
            {
                forAll(popBal_.sizeGroups(), i)
                {
                    const Foam::diameterModels::sizeGroup& fi =
                        popBal_.sizeGroups()[i];

                    scalar result(0.0);
                    scalar delta(0.0);

                    switch (abszissaType_)
                    {
                        case atVolume:
                        {
                            delta = popBal_.v()[i+1].value()
                              - popBal_.v()[i].value();

                            break;
                        }

                        case atDiameter:
                        {
                            const scalar& formFactor =
                                fi.VelocityGroup().formFactor().value();

                            delta =
                                pow
                                (
                                    popBal_.v()[i+1].value()
                                   /formFactor,
                                    1.0/3.0
                                )
                              - pow
                                (
                                    popBal_.v()[i].value()
                                   /formFactor,
                                    1.0/3.0
                                );

                            break;
                        }
                    }

                    switch (functionType_)
                    {
                        case ftNdf:
                        {
                            if (normalize_ == true)
                            {
                                result = N_[i]/delta/sumN_;
                            }
                            else
                            {
                                result = N_[i]/delta;
                            }

                            break;
                        }

                        case ftVdf:
                        {
                            if (normalize_ == true)
                            {
                                result = N_[i]*fi.x().value()/delta/sumV_;
                            }
                            else
                            {
                                result = N_[i]*fi.x().value()/delta;
                            }

                            break;
                        }

                        case ftNc:
                        {
                            if (normalize_ == true)
                            {
                                result = N_[i]/sumN_;
                            }
                            else
                            {
                                result = N_[i];
                            }

                            break;
                        }

                        default:
                        {
                            break;
                        }
                    }

                    file()<< tab << result;
                }
            }
        }
    }
    {
        file()<< endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
