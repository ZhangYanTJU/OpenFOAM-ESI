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

\*---------------------------------------------------------------------------*/

#include "fieldStatistics.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldStatistics, 0);
    addToRunTimeSelectionTable(functionObject, fieldStatistics, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::fieldStatistics::modeType
>
Foam::functionObjects::fieldStatistics::modeTypeNames_
({
    { modeType::MAG,  "magnitude" },
    { modeType::CMPT, "component" },
});

const Foam::Enum
<
    Foam::functionObjects::fieldStatistics::meanType
>
Foam::functionObjects::fieldStatistics::meanTypeNames_
({
    { meanType::ARITHMETIC, "arithmetic" },
    { meanType::VOLUMETRIC, "volumetric" },
});

const Foam::Enum
<
    Foam::functionObjects::fieldStatistics::calcType
>
Foam::functionObjects::fieldStatistics::calcTypeNames_
({
    // UNKNOWN is not enumerated
    { calcType::MIN, "min" },
    { calcType::MAX, "max" },
    { calcType::MEAN, "mean" },
    { calcType::VARIANCE, "variance" },
});


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Implementation
#include "fieldStatisticsImpl.cxx"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::functionObjects::fieldStatistics::statistic
Foam::functionObjects::fieldStatistics::createStatistic
(
    const word& statName,
    const modeType mode
)
{
    statistic stat;
    stat.name_ = statName;

    const auto statType(calcTypeNames_.lookup(statName, calcType::UNKNOWN));

    stat.calc = [this, statType, mode](variantInput input) -> variantOutput
    {
        return std::visit
        (
            [this, statType, mode](auto&& arg) -> variantOutput
            {
                using T = std::decay_t<decltype(arg)>;
                using value_type = typename T::value_type;

                switch (statType)
                {
                    case calcType::MIN :
                    {
                        if (mode == modeType::MAG)
                            return calcMin<scalar>(mag(arg));
                        else
                            return calcMin<value_type>(arg);
                    }
                    case calcType::MAX :
                    {
                        if (mode == modeType::MAG)
                            return calcMax<scalar>(mag(arg));
                        else
                            return calcMax<value_type>(arg);
                    }
                    case calcType::MEAN :
                    {
                        if (mode == modeType::MAG)
                            return calcMean<scalar>(mag(arg));
                        else
                            return calcMean<value_type>(arg);
                    }
                    case calcType::VARIANCE :
                    {
                        if (mode == modeType::MAG)
                            return calcVariance<scalar>(mag(arg));
                        else
                            return calcVariance<value_type>(arg);
                    }
                    default :
                    {
                        // Default case (for compiler)
                        return scalar(0);
                    }
                }
            },
            input
        );
    };

    return stat;
}


void Foam::functionObjects::fieldStatistics::writeFileHeader
(
    Ostream& os,
    const word& fieldName
)
{
    writeHeader(os, word("Field Statistics: " + fieldName));
    writeCommented(os, "Time");

    // Number of input statistics (i.e., statistics_) should be the same with
    // that of output statistics (i.e., results_). However, for consistency,
    // the output file columns are based on output statistics.
    const auto& result = results_(fieldName);

    for (const auto& iter : result.csorted())
    {
        const word& name = iter.key();
        writeTabbed(os, name);
    }
    os  << endl;
}


void Foam::functionObjects::fieldStatistics::writeStatData()
{
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        const auto& results = results_(fieldName);

        if (!results.size()) break;

        OFstream& file = *filePtrs_(fieldName);

        writeCurrentTime(file);

        for (const auto& iter : results.csorted())
        {
            const variantOutput& value = iter.val();

            std::visit
            (
                [&file](const auto& v)
                {
                    if constexpr
                    (
                        is_vectorspace_v<std::decay_t<decltype(v)>>
                    )
                    {
                        for (const auto& val : v) file<< token::TAB << val;
                    }
                    else
                    {
                        file<< token::TAB << v;
                    }
                },
                value
            );
        }
        file<< nl;
    }
}


void Foam::functionObjects::fieldStatistics::logStatData()
{
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        const auto& results = results_(fieldName);

        if (!results.size()) break;

        const word outputName
        (
            (mode_ == modeType::MAG)
          ? word("mag(" + fieldName + ")")
          : fieldName
        );

        Info<< nl << "    Field " << outputName << nl;

        for (const auto& iter : results.csorted())
        {
            const word& name = iter.key();
            const variantOutput& value = iter.val();

            Info<< "    " << name;
            std::visit
            (
                [](const auto& v)
                {
                    if constexpr
                    (
                        is_vectorspace_v<std::decay_t<decltype(v)>>
                    )
                    {
                        for (const auto& val : v) Info<< ' ' << val;
                    }
                    else
                    {
                        Info<< ' ' << v;
                    }
                },
                value
            );
            Info<< nl;
        }
    }
    Info<< endl;
}


void Foam::functionObjects::fieldStatistics::writeExtremaFileHeader
(
    Ostream& os,
    const word& fieldName
)
{
    writeHeader(os, word("Field Extrema Data: " + fieldName));
    writeCommented(os, "Time");
    writeTabbed(os, "min");
    writeTabbed(os, "min_procID");
    writeTabbed(os, "min_cellID");
    writeTabbed(os, "min_position");
    writeTabbed(os, "max");
    writeTabbed(os, "max_procID");
    writeTabbed(os, "max_cellID");
    writeTabbed(os, "max_position");
    os  << endl;
}


void Foam::functionObjects::fieldStatistics::writeExtremaData()
{
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        const auto& min = extremaResults_(fieldName).first();
        const auto& max = extremaResults_(fieldName).second();

        OFstream& file = *extremaFilePtrs_(fieldName);

        writeCurrentTime(file);

        file<< token::TAB;

        std::visit([&file](const auto& v){ file<< v; }, min.value_);

        file<< token::TAB
            << min.procID_ << token::TAB
            << min.cellID_ << token::TAB
            << min.position_ << token::TAB;

        std::visit([&file](const auto& v){ file<< v; }, max.value_);

        file<< token::TAB
            << max.procID_ << token::TAB
            << max.cellID_ << token::TAB
            << max.position_ << nl;
    }
}


void Foam::functionObjects::fieldStatistics::logExtremaData()
{
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        const auto& min = extremaResults_(fieldName).first();
        const auto& max = extremaResults_(fieldName).second();

        const word outputName
        (
            (mode_ == modeType::MAG)
          ? word("mag(" + fieldName + ")")
          : fieldName
        );

        std::visit
        (
            [outputName](const auto& v)
            {
                Info<< "    min(" << outputName << ") = " << v;
            },
            min.value_
        );

        Info<< " in cell " << min.cellID_
            << " at location " << min.position_
            << " on processor " << min.procID_;

        std::visit
        (
            [outputName](const auto& v)
            {
                Info<< nl << "    max(" << outputName << ") = " << v;
            },
            max.value_
        );

        Info<< " in cell " << max.cellID_
            << " at location " << max.position_
            << " on processor " << max.procID_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldStatistics::fieldStatistics
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    internal_(false),
    extrema_(false),
    mode_(modeType::MAG),
    mean_(meanType::ARITHMETIC),
    fieldSet_(mesh_)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldStatistics::read(const dictionary& dict)
{
    if (!(fvMeshFunctionObject::read(dict) && writeFile::read(dict)))
    {
        return false;
    }

    internal_ = dict.getOrDefault("internal", false);
    extrema_ = dict.getOrDefault("extrema", false);

    mode_ = modeTypeNames_.getOrDefault("mode", dict, modeType::MAG);
    mean_ = meanTypeNames_.getOrDefault("mean", dict, meanType::ARITHMETIC);

    // Reset and reprepare the input field names
    fieldSet_.clear();
    fieldSet_.read(dict);

    // Reset and reprepare the input statistics
    statistics_.clear();

    const wordHashSet statNames(dict.get<wordHashSet>("statistics"));
    // NOTE: could also filter out bad/unknown types here

    for (const word& statName : statNames.sortedToc())
    {
        statistics_.insert(statName, createStatistic(statName, mode_));
    }

    // Reset the output-statistics container
    results_.clear();

    // Reset the extrema-data container
    extremaResults_.clear();

    return true;
}


bool Foam::functionObjects::fieldStatistics::execute()
{
    fieldSet_.updateSelection();

    // Calculate the specified statistics for the specified fields
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        const bool ok =
        (
            calcStat<scalar>(fieldName)
         || calcStat<vector>(fieldName)
         || calcStat<sphericalTensor>(fieldName)
         || calcStat<symmTensor>(fieldName)
         || calcStat<tensor>(fieldName)
        );

        if (!ok)
        {
            WarningInFunction
                << "Unable to find field " << fieldName << endl;
        }
    }

    // Store the statistical results into the state containers
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        const auto& results = results_(fieldName);

        for (const auto& iter : results.csorted())
        {
            const word& name = iter.key();
            const variantOutput& value = iter.val();

            const word variableName(fieldName + "_" + name);

            std::visit
            (
                [this, variableName](const auto& v)
                {
                    this->setResult(variableName, v);
                },
                value
            );
        }

        if (extrema_)
        {
            const auto& min = extremaResults_(fieldName).first();
            std::visit
            (
                [this, fieldName](const auto& v)
                {
                    this->setResult(word(fieldName + "_min"), v);
                },
                min.value_
            );
            this->setResult(word(fieldName + "_min_procID"), min.procID_);
            this->setResult(word(fieldName + "_min_cellID"), min.cellID_);
            this->setResult(word(fieldName + "_min_position"), min.position_);

            const auto& max = extremaResults_(fieldName).second();
            std::visit
            (
                [this, fieldName](const auto& v)
                {
                    this->setResult(word(fieldName + "_max"), v);
                },
                max.value_
            );
            this->setResult(word(fieldName + "_max_procID"), max.procID_);
            this->setResult(word(fieldName + "_max_cellID"), max.cellID_);
            this->setResult(word(fieldName + "_max_position"), max.position_);
        }
    }

    return true;
}


bool Foam::functionObjects::fieldStatistics::write()
{
    Log << type() << ' ' << name() << " write:" << endl;

    // Create an output file per field
    if (writeToFile() && !writtenHeader_)
    {
        for (const word& fieldName : fieldSet_.selectionNames())
        {
            filePtrs_.set(fieldName, newFileAtStartTime(fieldName));

            OFstream& file = *filePtrs_(fieldName);

            writeFileHeader(file, fieldName);
        }

        if (extrema_)
        {
            for (const word& fieldName : fieldSet_.selectionNames())
            {
                extremaFilePtrs_.set
                (
                    fieldName,
                    newFileAtStartTime(word("fieldMinMax_" + fieldName))
                );

                OFstream& file = *extremaFilePtrs_(fieldName);
                writeExtremaFileHeader(file, fieldName);
            }
        }

        writtenHeader_ = true;
    }

    // Write the statistical results to the output file if requested
    if (writeToFile())
    {
        if (extrema_) writeExtremaData();

        writeStatData();
    }

    // Print the statistical results to the standard stream if requested
    if (log)
    {
        if (extrema_) logExtremaData();

        logStatData();
    }

    return true;
}


// ************************************************************************* //
