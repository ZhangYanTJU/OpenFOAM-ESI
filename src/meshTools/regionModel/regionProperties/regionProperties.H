/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2025 OpenCFD Ltd.
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

Class
    Foam::regionProperties

Description
    Simple class to hold region information for coupled region simulations.

    Gives per physics ('fluid', 'solid') the names of the regions. There
    is no assumption on this level that a region should only have a single
    set of physics.

    Uses the "regions" table from the constant/regionProperties file.
    For example,

    \verbatim
    regions
    (
        fluid   (air water)
        solid   (walls)
    );
    \endverbatim

SourceFiles
    regionProperties.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_regionProperties_H
#define Foam_regionProperties_H

#include "HashTable.H"
#include "fileName.H"
#include "wordList.H"
#include "IOobjectOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;

/*---------------------------------------------------------------------------*\
                       Class regionProperties Declaration
\*---------------------------------------------------------------------------*/

class regionProperties
:
    public HashTable<wordList>
{
public:

    // Constructors

        //- Construct from Time with specified read option
        //- (default: MUST_READ)
        explicit regionProperties
        (
            const Time& runTime,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ
        );

        //- Construct from Time and local (prefix) with specified read option
        //- (default: MUST_READ)
        regionProperties
        (
            const Time& runTime,
            const fileName& local,
            IOobjectOption::readOption rOpt = IOobjectOption::MUST_READ
        );


    //- Destructor
    ~regionProperties() = default;


    // Member Functions

        //- Total count of all region names.
        label count() const;

        //- The region names. Sorted by region type.
        wordList names() const;

        //- The region names in sorted order.
        wordList sortedNames() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
