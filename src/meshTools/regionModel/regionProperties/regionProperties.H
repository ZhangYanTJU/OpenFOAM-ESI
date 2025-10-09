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
    regionProperties.cxx
    regionModelProperties.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_regionProperties_H
#define Foam_regionProperties_H

#include "fileName.H"
#include "wordList.H"
#include "HashTable.H"
#include "IOobjectOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// General controls

namespace Foam
{
namespace regionModels
{

//- The enable/disable state for regionFaModel (default: true)
bool allowFaModels() noexcept;

//- Set enable/disable state for regionFaModel
//  \return previous state
bool allowFaModels(bool on) noexcept;

} // End namespace regionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class wordRes;

/*---------------------------------------------------------------------------*\
                       Class regionProperties Declaration
\*---------------------------------------------------------------------------*/

class regionProperties
:
    public HashTable<wordList>
{
public:

    // Generated Methods

        //- Default construct
        regionProperties() = default;

        //- Copy construct
        regionProperties(const regionProperties&) = default;

        //- Move construct
        regionProperties(regionProperties&&) = default;

        //- Copy assignment
        regionProperties& operator=(const regionProperties&) = default;

        //- Move assignment
        regionProperties& operator=(regionProperties&&) = default;


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


    // Static Member Functions

        //- The expected search path (eg, constant/regionProperties).
        //  Can be useful when generation error messages
        static fileName objectRelPath
        (
            const Time& runTime,
            const fileName& local = fileName::null
        );


    // Member Functions

        //- Total count of all region names.
        label count() const;

        //- The names of the region types (sorted)
        auto types() const { return HashTable<wordList>::sortedToc(); }

        //- The region names (sorted by region type)
        wordList names() const;

        //- The region names in sorted order.
        wordList sortedNames() const;

        //- The region names (sorted by region type) for all matches
        //  \returns an empty list for an empty matcher
        wordList names(const wordRes& matcher) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
