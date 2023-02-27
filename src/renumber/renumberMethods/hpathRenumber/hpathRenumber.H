/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 Alon Zameret, Noam Manaker Morag
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
    Foam::hpathRenumber

Description
    Renumbering based on Hamiltonian path.

SourceFiles
    hpathRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_hpathRenumber_H
#define Foam_hpathRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class hpathRenumber Declaration
\*---------------------------------------------------------------------------*/

class hpathRenumber
:
    public renumberMethod
{
    // Private Data

        //- Apply layer separation? (default: true)
        bool applyLayerSeparation_;


    // Private Classes

        // Helper class
        class hpathFinder;

public:

    // Generated Methods

        //- No copy construct
        hpathRenumber(const hpathRenumber&) = delete;

        //- No copy assignment
        void operator=(const hpathRenumber&) = delete;


    //- Runtime type information
    TypeName("hpath");


    // Constructors

        //- Construct given the renumber dictionary
        explicit hpathRenumber(const dictionary& dict);


    //- Destructor
    virtual ~hpathRenumber() = default;


    // Member Functions

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        virtual labelList renumber(const pointField&) const
        {
            NotImplemented;
            return labelList();
        }

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        //  Uses mesh for regIOobject
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField& cc
        ) const;

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        virtual labelList renumber
        (
            const CompactListList<label>& cellCells,
            const pointField& cellCentres
        ) const
        {
            NotImplemented;
            return labelList();
        }

        //- Return the order in which cells need to be visited
        //- (ie. from ordered back to original cell label).
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cellCentres
        ) const
        {
            NotImplemented;
            return labelList();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
