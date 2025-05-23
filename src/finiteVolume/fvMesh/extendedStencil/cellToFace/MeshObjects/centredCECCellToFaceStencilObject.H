/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::centredCECCellToFaceStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCECCellToFaceStencilObject_H
#define centredCECCellToFaceStencilObject_H

#include "extendedCentredCellToFaceStencil.H"
#include "CECCellToFaceStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredCECCellToFaceStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCECCellToFaceStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCECCellToFaceStencilObject
    >,
    public extendedCentredCellToFaceStencil
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            TopologicalMeshObject,
            centredCECCellToFaceStencilObject
        > MeshObject_type;

public:

    TypeName("centredCECCellToFaceStencil");


    // Constructors

        //- Construct from uncompacted face stencil
        explicit centredCECCellToFaceStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject_type(mesh),
            extendedCentredCellToFaceStencil(CECCellToFaceStencil(mesh))
        {
            if (extendedCellToFaceStencil::debug)
            {
                Info<< "Generated centred stencil " << type()
                    << nl << endl;
                writeStencilStats(Info, stencil(), map());
            }
        }


    //- Destructor
    virtual ~centredCECCellToFaceStencilObject() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
