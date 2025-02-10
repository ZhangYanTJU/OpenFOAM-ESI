/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2025 OpenCFD Ltd.
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

#include "simpleObjectRegistry.H"
#include "dictionary.H"
#include "ITstream.H"
#include "SpanStream.H"
#include "StringStream.H"
#include "int.H"
#include "floatScalar.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleObjectRegistry::setValues
(
    const dictionary& dict,
    bool verbose,
    bool dryrun
)
{
    // Report enables output, but respect DetailInfo state as well.
    verbose = (verbose && Foam::infoDetailLevel > 0);

    if (dryrun && !verbose)
    {
        return;
    }

    for (const entry& dEntry : dict)
    {
        const word& name = dEntry.keyword();

        simpleObjectRegistryEntry* objPtr = this->find(name);

        if (verbose)
        {
            if (objPtr)
            {
                Info<< "    " << dEntry << nl;
            }
            else
            {
                Info<< "    " << name << " (unregistered)" << nl;
            }
        }

        if (dryrun)
        {
            // Nothing else to do
        }
        else if (objPtr)
        {
            const List<simpleRegIOobject*>& objects = *objPtr;

            OCharStream os;
            ISpanStream is;

            if (dEntry.isDict())
            {
                os.rewind();
                os << dEntry.dict();

                is.reset(os.view());

                // Or alternatively?
                // ITstream is(dEntry.dict().tokens());

                for (simpleRegIOobject* obj : objects)
                {
                    is.rewind();
                    obj->readData(is);
                }
            }
            else  // dEntry.isStream()
            {
                for (simpleRegIOobject* obj : objects)
                {
                    obj->readData(dEntry.stream());
                }
            }
        }
    }
}


void Foam::simpleObjectRegistry::setNamedValue
(
    const std::string_view name,
    int val,
    bool verbose,
    bool dryrun
)
{
    // Respect DetailInfo state
    verbose = (verbose && Foam::infoDetailLevel > 0);

    if (dryrun && !verbose)
    {
        return;
    }

    token tok(static_cast<label>(val));

    // Handle name=value,
    // treating 'name=' like 'name' (ie, default value)

    const auto eq = name.find('=');
    std::string_view objName = name;
    std::string_view param;
    std::string key;

    if (eq != std::string_view::npos)
    {
        key = name.substr(0, eq);
        param = name.substr(eq+1);
        objName = key;
    }

    simpleObjectRegistryEntry* objPtr = this->find(objName.data());

    // Fail early
    if (!objPtr)
    {
        if (verbose)
        {
            Info<< objName.data() << " (unregistered)" << nl;
        }
        return;
    }


    if (!param.empty())
    {
        float fvalue(0);

        if (Foam::readInt(param.data(), val))
        {
            // Parses as int
            tok = static_cast<label>(val);
        }
        else if (Foam::readFloat(param.data(), fvalue))
        {
            // Parses as float
            tok = fvalue;
        }
        else
        {
            // Accept 'name=string' for named enums
            tok = Foam::string(param.data(), param.size());
        }
    }

    if (verbose)
    {
        Info<< objName.data() << '=' << tok << nl;
    }

    if (dryrun)
    {
        // Nothing else to do
    }
    else if (objPtr)
    {
        // The generic interface requires an Istream.
        ITstream is(tokenList(Foam::one{}, std::move(tok)));

        const List<simpleRegIOobject*>& objects = *objPtr;

        for (simpleRegIOobject* obj : objects)
        {
            is.rewind();
            obj->readData(is);
        }
    }
}


// ************************************************************************* //
