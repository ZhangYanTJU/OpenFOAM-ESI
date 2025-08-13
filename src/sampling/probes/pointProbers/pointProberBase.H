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
    Base class for sampling fields at specified internal and boundary locations.

SourceFiles
    pointProberBase.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointProberBase_H
#define Foam_pointProberBase_H

#include "fvMesh.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pointProberBase Declaration
\*---------------------------------------------------------------------------*/

class pointProberBase
:
    public pointField
{
protected:

    template<class T>
    struct isNotEqOp
    {
        void operator()(T& x, const T& y) const
        {
            const T unsetVal(-VGREAT*pTraits<T>::one);

            if (x != unsetVal)
            {
                // Keep x.

                // Note: should check for y != unsetVal but multiple sample cells
                // already handled in read().
            }
            else
            {
                // x is not set. y might be.
                x = y;
            }
        }
    };


    // Protected Data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Fixed locations (default: true)
        //  Note: set to false for moving mesh calculations where locations
        //        should move with the mesh
        bool fixedLocations_;

        //- Include probes that were not found (default: true)
        bool includeOutOfBounds_;

        //- Interpolation/sample scheme to obtain values at the points
        //  Note: only possible when fixedLocations_ is true
        word samplePointScheme_;


    // Calculated

        //- Cells to be probed (obtained from the locations)
        labelList elementList_;

        //- Faces to be probed
        labelList faceList_;

        //- Processor holding the cell or face (-1 if point not found
        //- on any processor)
        labelList processor_;

        //- Patch IDs on which the new probes are located
        labelList patchIDList_;

        //- Original probes location
        pointField oldPoints_;


    // Protected Member Functions

        //- Find cells and faces containing probes
        virtual void findElements(const fvMesh& mesh) = 0;


public:

    //- Runtime type information
    TypeName("pointProberBase");


    // Generated Methods

        //- No copy construct
        pointProberBase(const pointProberBase&) = delete;

        //- No copy assignment
        void operator=(const pointProberBase&) = delete;


    // Constructors

        //- Construct from Time and dictionary
        pointProberBase
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~pointProberBase() = default;


    // Member Functions

    // Access

        //- Return true if no probe locations
        bool empty() const { return pointField::empty(); }

        //- Return true if fixed locations
        bool fixedLocations() const { return fixedLocations_; }

        //- Return true if include out of bounds probes
        bool includeOutOfBounds() const { return includeOutOfBounds_; }

        //- Return the interpolation scheme to obtain values at the points
        //  Note: only possible when fixedLocations_ is true
        const word& samplePointScheme() const { return samplePointScheme_; }

        //- Return const reference to the probe locations
        const pointField& probeLocations() const { return *this; }

        //- Return the location of probe i
        const point& probe(const label i) const { return operator[](i); }

        //- Cells to be probed (obtained from the locations)
        const labelList& elements() const { return elementList_; }

        //- Return const reference to the faces to be probed
        const labelList& faces() const { return faceList_; }

        //- Return const reference to the processor list
        const labelList& processors() const { return processor_; }

        //- Return const reference to the patch ID list
        const labelList& patchIDList() const noexcept { return patchIDList_; }

        //- Return const reference to the original probe locations
        const pointField& oldPoints() const noexcept { return oldPoints_; }


    // I-O

        //- Read the settings dictionary
        virtual bool read(const dictionary&);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
