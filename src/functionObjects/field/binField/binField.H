/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::functionObjects::binField

Group
    grpFieldFunctionObjects

Description
Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    \endverbatim

    where the entries mean:
    \table
      Property   | Description               | Type     | Rqd | Dflt
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

See also
  - Foam::functionObject
  - Foam::functionObjects::fvMeshFunctionObject
  - Foam::functionObjects::writeFile

SourceFiles
    binField.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_binField_H
#define functionObjects_binField_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class binModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class binField Declaration
\*---------------------------------------------------------------------------*/

class binField
:
    public fvMeshFunctionObject
{
protected:

    // Protected Data

        //- Runtime-selectable bin model
        autoPtr<binModel> binModelPtr_;


public:

    //- Runtime type information
    TypeName("binField");


    // Constructors

        //- No default construct
        binField() = delete;

        //- Construct from Time and dictionary
        binField
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool readFields = true
        );

        //- Construct from objectRegistry and dictionary
        binField
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool readFields = true
        );

        //- No copy construct
        binField(const binField&) = delete;

        //- No copy assignment
        void operator=(const binField&) = delete;


    //- Destructor
    virtual ~binField() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);

        //- Execute the function object
        virtual bool execute();

        //- Write to data files/fields and to streams
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
