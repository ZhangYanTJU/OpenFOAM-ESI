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

Class
    Foam::zoneSubSet

Description
    Intermediate tool for handling \c cellZones
    for function objects (e.g. \c momentumError)
    wherein the sub-mesh option is available.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    <functionObject>
    {
        // Mandatory and optional entries
        ...

        // Inherited entries
        cellZones         ( <cellZone1> <cellZone2> ... <cellZoneN> );
        nLayers           <label>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type  | Reqd | Deflt
      cellZones | Names of cellZones                    | wordRes | no   | -
      nLayers   | Number of cell layers around cellZones  | label | no   | 0
    \endtable

See also
  - Foam::functionObjects::momentumError
  - Foam::functionObjects::div

SourceFiles
    zoneSubSet.C
    zoneSubSetTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef zoneSubSet_H
#define zoneSubSet_H

#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class zoneSubSet Declaration
\*---------------------------------------------------------------------------*/

class zoneSubSet
{
    // Private Data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Sub-set mesh
        autoPtr<fvMeshSubset> subSetMeshPtr_;

        //- Region (zone/surface) names
        wordRes zoneNames_;

        //- Number of layer around zones
        label nLayers_;


    // Private Members Functions

        //- Initialise sub-mesh environment
        void init(const dictionary& dict);


        //- No copy construct
        zoneSubSet(const zoneSubSet&) = delete;

        //- No copy assignment
        void operator=(const zoneSubSet&) = delete;


public:

    //- Runtime type information
    TypeName("zoneSubSet");


    // Constructors

        //- Construct from Time and dictionary
        zoneSubSet
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~zoneSubSet() = default;


    // Member Functions

        //- Map from the sub-mesh to original cell zones
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> mapToZone
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );


        // Access

            //- Return const access to the subSetMesh
            const fvMeshSubset& subSetMesh() const
            {
                return subSetMeshPtr_();
            }

            //- Return non-const access to the subSetMesh
            fvMeshSubset& subSetMesh()
            {
                return subSetMeshPtr_.ref();
            }

            //- Return const access to the zone names
            const wordRes& zoneNames() const noexcept
            {
                return zoneNames_;
            }

            //- Return number of layers around cell zones
            label nLayers() noexcept
            {
                return nLayers_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zoneSubSetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
