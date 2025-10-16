/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

Namespace
    Foam::cornerDetectionModels

Description
    A namespace for various \c cornerDetection model implementations.

Class
    Foam::cornerDetectionModel

Description
    A base class for \c cornerDetection models.

SourceFiles
    cornerDetectionModel.C
    cornerDetectionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cornerDetectionModel_H
#define Foam_cornerDetectionModel_H

#include "liquidFilmBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                    Class cornerDetectionModel Declaration
\*---------------------------------------------------------------------------*/

class cornerDetectionModel
{
    // Private Data

        //- Const reference to the finite-area mesh
        const faMesh& mesh_;

        //- Const reference to the film
        const regionModels::areaSurfaceFilmModels::liquidFilmBase& film_;

        //- Const reference to the dictionary
        const dictionary& dict_;

        //- Bitset for corner faces, true for corner faces
        bitSet cornerFaces_;

        //- List of corner angles [rad]
        scalarList cornerAngles_;


protected:

    // Protected Member Functions

        //- Return the dihedral angle between two neighbour faces
        //  Robust 2*atan form (handles parallel normals better than acos)
        //  \param n0 First normal vector
        //  \param n1 Second normal vector
        //  \return The dihedral angle [rad]
        scalar dihedralAngle(const vector& n0, const vector& n1) const;


public:

    //- Runtime type information
    TypeName("cornerDetectionModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cornerDetectionModel,
            dictionary,
            (
                const faMesh& mesh,
                const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
                const dictionary& dict
            ),
            (mesh, film, dict)
        );


    // Selectors

        //- Return a reference to the selected cornerDetection model
        static autoPtr<cornerDetectionModel> New
        (
            const faMesh& mesh,
            const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
            const dictionary& dict
        );


    // Constructors

        //- Construct from components
        cornerDetectionModel
        (
            const faMesh& mesh,
            const regionModels::areaSurfaceFilmModels::liquidFilmBase& film,
            const dictionary& dict
        );


    // Generated Methods

        //- No copy construct
        cornerDetectionModel(const cornerDetectionModel&) = delete;

        //- No copy assignment
        void operator=(const cornerDetectionModel&) = delete;


    //- Destructor
    virtual ~cornerDetectionModel();


    // Member Functions

    // Getters

        //- Return const reference to the finite-area mesh
        const faMesh& mesh() const noexcept { return mesh_; }

        //- Return const reference to the film
        const regionModels::areaSurfaceFilmModels::liquidFilmBase&
            film() const noexcept { return film_; }

        //- Return const reference to the dictionary
        const dictionary& dict() const noexcept { return dict_; }

        //- Return const reference to the corner faces bitset
        const bitSet& getCornerFaces() const noexcept { return cornerFaces_; }

        //- Return const reference to the corner angles list [rad]
        const scalarList& getCornerAngles() const noexcept
        {
            return cornerAngles_;
        }


    // Setters

        //- Set the corner faces bitset
        void setCornerFaces(const bitSet& cornerFaces)
        {
            cornerFaces_ = cornerFaces;
        }

        //- Set the corner angles list [rad]
        void setCornerAngles(const scalarList& cornerAngles)
        {
            cornerAngles_ = cornerAngles;
        }


    // Evaluation

        //- Detect and store the corner faces
        virtual bool detectCorners() = 0;


    // I-O

        //- Read the model dictionary
        virtual bool read(const dictionary&) { return true; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
