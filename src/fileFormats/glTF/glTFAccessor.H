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
    Foam::glTF::accessor

Description
    glTF accessor

Note
    Implements the glTF v2 specification

SourceFiles
    glTFAccessor.C
    glTFAccessorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef gltf_accessor_H
#define gltf_accessor_H

#include "glTFBase.H"
#include "Field.H"
#include "label.H"
#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace glTF
{
    class accessor;
}

Ostream& operator<<(Ostream& os, const glTF::accessor& acc);

namespace glTF
{

/*---------------------------------------------------------------------------*\
                       Class glTF::accessor Declaration
\*---------------------------------------------------------------------------*/

class accessor
:
    public base
{
protected:

    // Protected Data

        //- Buffer view index
        label bufferViewId_;

        //- Byte offset
        label byteOffset_;

        //- Component type
        label componentType_;

        //- Data size
        label count_;

        //- Type
        string type_;

        //- Max value.  Note: stored as a string for convenience
        string max_;

        //- Min value.  Note: stored as a string for convenience
        string min_;

        //- Flag to indicate whether min and max values are available
        bool minMax_;


public:

    // Constructors

        //- Null constructor
        accessor();

        //- Construct from name
        accessor(const word& name);


    //- Destructor
    ~accessor() = default;


    // Public Member Functions

        //- Return the buffer view index
        label& bufferViewId();

        //- Return the byte offset
        label& byteOffset();

        //- Return the component type
        label& componentType();

        //- Return the data size
        label& count();

        //- Return the type
        string& type();

        //- Helper function to return the glTF value type that matches the
        //- OpenFOAM primitive type
        template<class Type>
        string getValueType() const;

        //- Set the accessor
        template<class Type>
        void set(const Field<Type>& fld, bool minMax = true);

        //- Write
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace glTF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "glTFAccessorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //