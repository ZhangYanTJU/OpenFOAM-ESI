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
    Foam::multiWorldConnections

Description
    Centralized handling of multi-world MPI connections.

SourceFiles
    multiWorldConnectionsObject.C

\*---------------------------------------------------------------------------*/

#ifndef meshObjects_multiWorldConnections_H
#define meshObjects_multiWorldConnections_H

#include "MeshObject.H"
#include "EdgeMap.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiWorldConnections Declaration
\*---------------------------------------------------------------------------*/

class multiWorldConnections
:
    public MeshObject
    <
        Time,
        TopologicalMeshObject,
        multiWorldConnections
    >
{
    // Private Typedefs

        typedef MeshObject
        <
            Time,
            TopologicalMeshObject,
            multiWorldConnections
        > MeshObjectType;


    // Private Data

        //- Table of world/world connections to communicator label
        EdgeMap<label> table_;


    // Private Member Functions

        //- Pairing from myWorld to other world by ID
        //  \return invalid on any error or if non-parallel
        static edge worldPair(const label otherWorld);

        //- Pairing from myWorld to other world by NAME.
        //  \return invalid on any error or if non-parallel
        static edge worldPair(const word& otherWorld);

        //- Allocate a communicator between myWorld and other world
        static label communicator(const edge& worlds);


        //- Define a connection from myWorld to other world by NAME
        bool addConnectionByName(const word& otherWorld);

        //- Get communicator for myWorld to other world connection by NAME.
        label getCommunicatorByName(const word& otherWorld) const;


public:

    //- Run-time type information
    TypeName("multiWorld");


    // Constructors

        //- Construct
        explicit multiWorldConnections(const Time& runTime);


    // Selectors

        //- Access mesh object
        static const multiWorldConnections& New(const Time& runTime);


    //- Destructor
    ~multiWorldConnections();


    // Member Functions

        //- Create all communicators.
        //- Low-level, not normally called directly
        void createComms();

        //- Define a connection from myWorld to other world by NAME
        bool addConnection(const word& otherWorld);

        //- Define a connection from myWorld to other world by ID
        bool addConnectionById(const label otherWorld);

        //- Get communicator for myWorld to other world connection by NAME.
        //  Uses cached value, or creates a new communicator
        label getCommunicator(const word& otherWorld) const;

        //- Get communicator for myWorld to other world connection by ID.
        //  Uses cached value, or creates a new communicator
        label getCommunicatorById(const label otherWorld) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
