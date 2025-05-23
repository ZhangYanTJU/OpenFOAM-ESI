/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "removeEntry.H"
#include "dictionary.H"
#include "wordRes.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        removeEntry,
        execute,
        dictionaryIstream,
        remove
    );
} // End namespace functionEntries
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::removeEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    const wordRes patterns(functionEntry::readStringList<wordRe>(is));

    for (const wordRe& key : patterns)
    {
        if (key.isLiteral() && key.contains('/'))
        {
            // Remove scoped keyword, or keyword in the local scope
            auto finder(parentDict.searchScoped(key, keyType::LITERAL));

            if (finder.good())
            {
                finder.context().remove(finder.ptr()->keyword());
            }
        }
        else
        {
            // Remove by pattern
            const wordList dictKeys = parentDict.toc();
            const labelList indices = wordRes::matching(key, dictKeys);

            for (const auto idx : indices)
            {
                parentDict.remove(dictKeys[idx]);
            }
        }
    }

    return true;
}


// ************************************************************************* //
