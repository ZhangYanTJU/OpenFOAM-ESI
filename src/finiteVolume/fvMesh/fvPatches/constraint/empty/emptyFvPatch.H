/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::emptyFvPatch

Description
    A patch which will not exist in the fvMesh. Typical example is a front and
    back plane of a 2-D geometry

SourceFiles
    emptyFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_emptyFvPatch_H
#define Foam_emptyFvPatch_H

#include "fvPatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class emptyFvPatch Declaration
\*---------------------------------------------------------------------------*/

class emptyFvPatch
:
    public fvPatch
{
    // Private Data

        //- Empty face-cell addressing
        const labelList::subList faceCells_;


public:

    //- Runtime type information
    TypeName(emptyPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        emptyFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm);


    // Member Functions

        //- Empty patch has zero size
        virtual label size() const { return 0; }

        //- Return faceCells
        virtual const labelUList& faceCells() const { return faceCells_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
