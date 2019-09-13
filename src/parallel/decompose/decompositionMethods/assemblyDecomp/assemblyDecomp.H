/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Foam::assemblyDecomp

Description
    Domain decomposition that conserves mapped patches locally in
    processors. This is necessary in cases where the coupled option
    for energy is used in order to make internal all faces
    of the mapped patch.

    This method uses a run-time selectable method for the designated main
    region. This decomposition will write the decomposition list for each
    of the remaining region which need to use 'manual' decomposition.

    \verbatim
        numberOfSubdomains  4;
        method  assembly;
        coeffs
        {
            numberOfSubdomains  4;
            method              scotch;
            dataFile            "cellDecompose.dat";
        }
    \endverbatim

SourceFiles
    assemblyDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef assemblyDecomp_H
#define assemblyDecomp_H

#include "decompositionMethod.H"
#include "mappedPatchBase.H"
#include "fvMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class assemblyDecomp Declaration
\*---------------------------------------------------------------------------*/

class assemblyDecomp
:
    public decompositionMethod
{
    // Private Data

        //- Method dictionary
        dictionary methodDict_;

        //- Method
        autoPtr<decompositionMethod> method_;

        //- Data file for regions
        fileName dataFile_;


    // Private Member Functions

        //- No copy construct
        void operator=(const assemblyDecomp&) = delete;

        //- No copy assignment
        assemblyDecomp(const assemblyDecomp&) = delete;


public:

    //- Runtime type information
    TypeName("assembly");


    // Constructors

        //- Construct given the decomposition dictionary
        assemblyDecomp(const dictionary& decompDict);


    //- Destructor
    virtual ~assemblyDecomp() = default;


    // Member Functions

        //- This decompose does not care about proc boundaries.
        //- Is all up to the user.
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        ) const;

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
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
