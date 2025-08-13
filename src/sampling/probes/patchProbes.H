/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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
    Foam::patchProbes

Description
    Probes the specified points on specified patches. The points get snapped
    onto the nearest point on the nearest face of the specified patch, and the
    sampling is actioned on the snapped locations.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    patchProbes
    {
        // Mandatory entries
        type            patchProbes;
        libs            (sampling);

        fields          (<wordRes>);
        probeLocations  (<vectorList>);
        patches         (<wordRes>);    // or patch  <word>;

        // Optional entries
        verbose         <bool>;
        sampleOnExecute <bool>;
        fixedLocations  <bool>;
        includeOutOfBounds <bool>;
        interpolationScheme <word>;

        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: patchProbes             | word | yes  | -
      libs         | Library name: sampling             | word | yes  | -
      fields       | Names of the fields to be probed   | wordRes | yes | -
      probeLocations | Locations of the probes          | vectorField | yes | -
      patches      | Patches to sample (wildcards allowed) | wordRes | yes | -
      verbose      | Enable/disable verbose output         | bool | no | false
      sampleOnExecute | Sample on execution and store results | bool | no <!--
                   -->                                               | false
      fixedLocations | Do not recalculate cells if mesh moves | bool | no | true
      includeOutOfBounds | Include out-of-bounds locations | bool | no | true
      interpolationScheme | Scheme to obtain values at the points | word <!--
                   -->                                            | no | cell
    \endtable

    The inherited entries are elaborated in:
     - \link fvMeshFunctionObject.H \endlink
     - \link ProbesBase.H \endlink
     - \link patchPointProber.H \endlink
     - \link pointProberBase.H \endlink

Note
  - The \c includeOutOfBounds filters out points that haven't been found.
    Default is to include them (with value \c -VGREAT).

SourceFiles
    patchProbes.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_patchProbes_H
#define Foam_patchProbes_H

#include "ProbesBase.H"
#include "patchPointProber.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class patchProbes Declaration
\*---------------------------------------------------------------------------*/

class patchProbes
:
    public ProbesBase<patchPointProber>
{
    // Private Data

        //- Use simpler synonym for the base type
        using Base = ProbesBase<patchPointProber>;


public:

    //- Runtime type information
    TypeName("patchProbes");


    // Constructors

        //- Construct from name, Time and dictionary
        patchProbes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );


    //- Destructor
    virtual ~patchProbes() = default;


    // Member Functions

        //- Bring Base::prober into this class's public scope.
        using Base::prober;

        //- Read the settings from the dictionary
        virtual bool read(const dictionary&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
