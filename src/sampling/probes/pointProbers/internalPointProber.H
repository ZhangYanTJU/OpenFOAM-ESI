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

Class
    Foam::internalPointProber

Description
    A utility class for probing field values at specified point locations
    within an \c fvMesh.

    The \c internalPointProber stores a list of 3D point coordinates and
    determines the corresponding mesh elements (cells or faces) that contain
    these points. It provides methods to sample volume or surface fields at
    the stored locations, with support for fixed or mesh-moving point
    coordinates.

    Features include:
      - Reading probe locations and settings from a dictionary
      - Support for fixed or moving locations (for dynamic mesh cases)
      - Optional inclusion of points that lie outside of the mesh domain
      - Selection of interpolation/sampling schemes for fixed locations
      - Sampling of volume and surface fields by name or by direct reference
      - Automatic update of element mapping when the mesh changes or moves

SourceFiles
    internalPointProber.C
    internalPointProberTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_internalPointProber_H
#define Foam_internalPointProber_H

#include "pointProberBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class internalPointProber Declaration
\*---------------------------------------------------------------------------*/

class internalPointProber
:
    public pointProberBase
{
protected:

    // Protected Member Functions

        //- Find cells and faces containing probes
        virtual void findElements(const fvMesh& mesh);


public:

    //- Runtime type information
    TypeName("internalPointProber");


    // Constructors

        //- Construct from Time and dictionary
        internalPointProber
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~internalPointProber() = default;


    // Member Functions

    // Sampling

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample(const VolumeField<Type>&) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sample(const SurfaceField<Type>&) const;

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample(const word& fieldName) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sampleSurfaceField(const word& fieldName) const;


    // I-O

        //- Read the settings dictionary
        virtual bool read(const dictionary&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "internalPointProberTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
