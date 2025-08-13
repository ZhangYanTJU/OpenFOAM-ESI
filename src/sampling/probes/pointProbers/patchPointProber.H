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
    Foam::patchPointProber

Description
    Utility class for probing specified points on user-selected boundary
    patches. The input points are projected onto the nearest point of the
    nearest face on the specified patch, ensuring sampling occurs at valid
    patch locations.

    The patchPointProber enables sampling of both volume and surface fields
    at these snapped locations. Patch selection is controlled via patch names or
    indices, and the class provides runtime selection and dictionary-driven
    configuration.

    Typical usage involves specifying patch names and probe locations in a
    dictionary, after which the class manages the mapping and sampling
    operations.

SourceFiles
    patchPointProber.C
    patchPointProberTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_patchPointProber_H
#define Foam_patchPointProber_H

#include "pointProberBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class patchPointProber Declaration
\*---------------------------------------------------------------------------*/

class patchPointProber
:
    public pointProberBase
{
protected:

    // Protected Data

        //- Names of the patches to sample
        wordRes patchNames_;

        //- Index of the patches to sample
        labelList patchIDs_;


    // Protected Member Functions

        //- Find elements containing patchPointProber
        virtual void findElements(const fvMesh& mesh);


public:

    //- Runtime type information
    TypeName("patchPointProber");


    // Constructors

        //- Construct from Time and dictionary
        patchPointProber
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchPointProber() = default;


    // Member Functions

    // Access

        //- Return the index of the patches to sample
        const labelList& patchIDs() const noexcept { return patchIDs_; }


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
    #include "patchPointProberTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
