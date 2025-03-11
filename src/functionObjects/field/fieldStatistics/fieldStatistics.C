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
    { modeType::mdMag,  "magnitude" },
    { modeType::mdCmpt, "component" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::functionObjects::fieldStatistics::statistic
Foam::functionObjects::fieldStatistics::createStatistic
(
    const word& statName,
    const modeType& mode
)
{
    statistic m;
    m.name_ = statName;

    m.calc = [this, statName, mode](variantInput input) -> variantOutput
    {
        return std::visit
        (
            [this, statName, mode](auto&& arg) -> variantOutput
            {
                using T = std::decay_t<decltype(arg)>;
                using BaseType = typename T::value_type;

                if (statName == "min")
                {
                    if (mode == mdMag) return mag(calcMin<BaseType>(arg));
                    else return calcMin<BaseType>(arg);
                }
                else if (statName == "max")
                {
                    if (mode == mdMag) return mag(calcMax<BaseType>(arg));
                    else return calcMax<BaseType>(arg);
                }
                else if (statName == "mean")
                {
                    if (mode == mdMag) return mag(calcMean<BaseType>(arg));
                    else return calcMean<BaseType>(arg);
                }
                else if (statName == "variance")
                {
                    if (mode == mdMag) return mag(calcVariance<BaseType>(arg));
                    else return calcVariance<BaseType>(arg);
                }
                else return scalar{};  // Default case (for compiler)
            },
            input
        );
    };

    return m;
}


void Foam::functionObjects::fieldStatistics::writeFileHeader
(
    Ostream& os,
    const word& fieldName
)
{
    writeHeader(os, word("Field Statistics: " + fieldName));
    writeCommented(os, "Time");

    // Number of input statistics (i.e., statistics_) could be lesser than that
    // of output statistics (i.e., results_; e.g., cell index of min value).
    // Therefore, the output file columns are based on output statistics.
    const auto& result = results_(fieldName);

    forAllConstIters(result, iter)
    {
        const word& name = iter.key();
        writeTabbed(os, name);
    }
    os  << endl;

    writtenHeader_ = true;
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

    mode_ = modeTypeNames_.getOrDefault("mode", dict, modeType::mdMag);

    // Reset and reprepare the input field names
    fieldSet_.clear();
    fieldSet_.read(dict);
    fieldSet_.updateSelection();

    // Reset and reprepare the input statistics
    statistics_.clear();
    const wordHashSet stats(dict.get<wordHashSet>("statistics"));
    for (const word& m : stats)
    {
        statistics_.insert(m, createStatistic(m, mode_));
    }

    // Reset and reprepare the output statistics
    results_.clear();

    return true;
}


bool Foam::functionObjects::fieldStatistics::execute()
{
    // Calculate the specified statistics for the specified fields
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        const bool ok =
        (
            calcStat<scalar>(fieldName)
         || calcStat<vector>(fieldName)
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

        forAllConstIters(results, iter)
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
    }

    return true;
}


bool Foam::functionObjects::fieldStatistics::write()
{
    // Create an output file per field
    if (writeToFile() && !writtenHeader_)
    {
        for (const word& fieldName : fieldSet_.selectionNames())
        {
            filePtrs_.set(fieldName, newFileAtStartTime(fieldName));

            OFstream& file = *filePtrs_(fieldName);

            writeFileHeader(file, fieldName);
        }
    }

    // Write the statistical results to the output file if requested
    if (writeToFile())
    {
        for (const word& fieldName : fieldSet_.selectionNames())
        {
            const auto& results = results_(fieldName);

            OFstream& file = *filePtrs_(fieldName);

            writeCurrentTime(file);

            forAllConstIters(results, iter)
            {
                const variantOutput& value = iter.val();

                std::visit
                (
                    [&file](const auto& v)
                    {
                        if constexpr (std::is_same_v<std::decay_t<decltype(v)>, scalar>)
                        {
                            file<< token::TAB << v;
                        }
                        else
                        {
                            for (const auto& val : v) file<< token::TAB << val;
                        }
                    },
                    value
                );
            }
            file<< nl;
        }
    }

    // Print the statistical results to the standard stream if requested
    if (log)
    {
        for (const word& fieldName : fieldSet_.selectionNames())
        {
            const auto& results = results_(fieldName);

            forAllConstIters(results, iter)
            {
                const word& name = iter.key();
                const variantOutput& value = iter.val();

                std::visit
                (
                    [name](const auto& v)
                    {
                        if constexpr (std::is_same_v<std::decay_t<decltype(v)>, scalar>)
                        {
                            Info<< name << " " << v;
                        }
                        else
                        {
                            Info<< name << " ";
                            for (const auto& val : v) Info<< val << " ";
                        }
                        Info<< nl;
                    },
                    value
                );
            }
        }
    }

    return true;
}


// ************************************************************************* //
