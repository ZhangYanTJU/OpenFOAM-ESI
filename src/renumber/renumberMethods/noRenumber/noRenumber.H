/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::noRenumber

Description
    A dummy renumber method, selected as \c none.

    Method coefficients: \a none

SourceFiles
    noRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_noRenumber_H
#define Foam_noRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class noRenumber Declaration
\*---------------------------------------------------------------------------*/

class noRenumber
:
    public renumberMethod
{
public:

    //- Runtime type information
    TypeNameNoDebug("none");


    // Constructors

        //- Default construct
        noRenumber();

        //- Construct given the renumber dictionary (unused)
        explicit noRenumber(const dictionary& dict);


    //- Destructor
    virtual ~noRenumber() = default;


    // Member Functions

        //- The renumbering method does not require a polyMesh
        virtual bool needs_mesh() const { return false; }


        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        virtual labelList renumber(const pointField&) const;

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        virtual labelList renumber
        (
            //! Mesh number of cells
            const polyMesh& mesh,
            //! \em ignored
            const pointField& cellCentres
        ) const;

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        virtual labelList renumber
        (
            //! Mesh connectivity for number of cells
            const CompactListList<label>& cellCells,
            //! \em ignored
            const pointField& cellCentres
        ) const;

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        virtual labelList renumber
        (
            //! Mesh connectivity for number of cells
            const labelListList& cellCells,
            //! \em ignored
            const pointField& cellCentres
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
