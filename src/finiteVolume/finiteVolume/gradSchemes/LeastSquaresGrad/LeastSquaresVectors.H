/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::fv::LeastSquaresVectors

Description
    Least-squares gradient scheme vectors

See also
    Foam::fv::LeastSquaresGrad

SourceFiles
    LeastSquaresVectors.C

\*---------------------------------------------------------------------------*/

#ifndef LeastSquaresVectors_H
#define LeastSquaresVectors_H

#include "extendedCentredCellToCellStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class LeastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

template<class Stencil>
class LeastSquaresVectors
:
    public MeshObject<fvMesh, MoveableMeshObject, LeastSquaresVectors<Stencil>>
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            MoveableMeshObject,
            LeastSquaresVectors<Stencil>
        > MeshObject_type;


    // Private Data

        //- Least-squares gradient vectors
        List<List<vector>> vectors_;


    // Private Member Functions

        //- Calculate Least-squares gradient vectors
        void calcLeastSquaresVectors();


public:

    // Declare name of the class and its debug switch
    TypeName("LeastSquaresVectors");


    // Constructors

        //- Construct given an fvMesh and the minimum determinant criterion
        explicit LeastSquaresVectors(const fvMesh& mesh);


    //- Destructor
    virtual ~LeastSquaresVectors();


    // Member Functions

        //- Return const reference to the stencil
        const extendedCentredCellToCellStencil& stencil() const
        {
            return Stencil::New(this->mesh_);
        }

        //- Return const reference to the least square vectors
        const List<List<vector>>& vectors() const
        {
            return vectors_;
        }

        //- Update the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LeastSquaresVectors.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
