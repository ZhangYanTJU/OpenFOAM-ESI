/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::centredCPCCellToCellStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCPCCellToCellStencilObject_H
#define centredCPCCellToCellStencilObject_H

#include "extendedCentredCellToCellStencil.H"
#include "CPCCellToCellStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredCPCCellToCellStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCPCCellToCellStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCPCCellToCellStencilObject
    >,
    public extendedCentredCellToCellStencil
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            TopologicalMeshObject,
            centredCPCCellToCellStencilObject
        > MeshObject_type;

public:

    TypeName("centredCPCCellToCellStencil");


    // Constructors

        //- Construct from uncompacted cell stencil
        explicit centredCPCCellToCellStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject_type(mesh),
            extendedCentredCellToCellStencil(CPCCellToCellStencil(mesh))
        {}


    //- Destructor
    virtual ~centredCPCCellToCellStencilObject() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
