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
    Foam::fv::patchCellsSource

Group
    grpFvOptionsSources

Description
    Source defined by a boundary condition applied to cell next to patches.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    patchCellsSource1
    {
        // Mandatory entries (unmodifiable)
        type           patchCellsSource;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd | Dflt
      type      | Type name: patchCellsSource       | word   | yes   | -

    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link option.H \endlink


SourceFiles
    patchCellsSource.C
    patchCellsSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef patchCellsSource_H
#define patchCellsSource_H

#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class patchCellsSource Declaration
\*---------------------------------------------------------------------------*/

class patchCellsSource
:
    public fv::option
{

private:

    // Private Data

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Name of operand velocity field
        word UName_;

        //- Name of operand energy field
        word hName_;

        //- Name of operand species field
        word speciesName_;

public:

    //- Runtime type information
    TypeName("patchCellsSource");


    // Constructors

        //- Construct from explicit source name and mesh
        patchCellsSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        patchCellsSource
        (
            const patchCellsSource&
        ) = delete;

        //- No copy assignment
        void operator=(const patchCellsSource&) = delete;


    //- Destructor
    ~patchCellsSource() = default;


    // Member Functions


        //- Add explicit contribution to compressible enthalpy equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
