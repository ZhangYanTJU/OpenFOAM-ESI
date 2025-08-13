/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2025 OpenCFD Ltd.
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
    Foam::ProbesBase

Description
    Base class for sampling fields at specified locations and writing to file.
    The locations are specified and determined in the derived class. The
    sampling is done using the specified point prober class.

SourceFiles
    ProbesBase.C
    ProbesBaseTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ProbesBase_H
#define Foam_ProbesBase_H

#include "fvMeshFunctionObject.H"
#include "polyMesh.H"
#include "HashPtrTable.H"
#include "OFstream.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "pointProberBase.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ProbesBase Declaration
\*---------------------------------------------------------------------------*/

template<class Prober>
class ProbesBase
:
    public functionObjects::fvMeshFunctionObject
{
protected:

    // Protected Data

        //- The specified point prober
        Prober prober_;


    // Protected Classes

        //- Grouping of field names by GeometricField type
        template<class GeoField>
        struct fieldGroup : public DynamicList<word> {};


    // Data Types

        //- Local control for sampling actions
        enum sampleActionType : unsigned
        {
            ACTION_NONE  = 0,
            ACTION_WRITE = 0x1,
            ACTION_STORE = 0x2,
            ACTION_ALL = 0xF
        };


    // Protected Data

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        //- Output verbosity
        bool verbose_;

        //- Perform sample actions on execute as well
        bool onExecute_;

        //- Requested names of fields to probe
        wordRes fieldSelection_;


    // Calculated

        //- Current list of field names selected for sampling
        DynamicList<word> selectedFieldNames_;

        //- Categorized scalar/vector/tensor volume fields
        fieldGroup<volScalarField> scalarFields_;
        fieldGroup<volVectorField> vectorFields_;
        fieldGroup<volSphericalTensorField> sphericalTensorFields_;
        fieldGroup<volSymmTensorField> symmTensorFields_;
        fieldGroup<volTensorField> tensorFields_;

        //- Categorized scalar/vector/tensor surface fields
        fieldGroup<surfaceScalarField> surfaceScalarFields_;
        fieldGroup<surfaceVectorField> surfaceVectorFields_;
        fieldGroup<surfaceSphericalTensorField> surfaceSphericalTensorFields_;
        fieldGroup<surfaceSymmTensorField> surfaceSymmTensorFields_;
        fieldGroup<surfaceTensorField> surfaceTensorFields_;

        //- Current open files (non-empty on master only)
        HashPtrTable<OFstream> probeFilePtrs_;


    // Protected Member Functions

        //- Classify field types, close/open file streams
        //  \return number of fields to sample
        label prepare(unsigned request);

        //- Get from registry or load from disk
        template<class GeoField>
        tmp<GeoField> getOrLoadField(const word& fieldName) const;

        //- Store results: min/max/average/size
        template<class Type>
        void storeResults(const word& fieldName, const Field<Type>& values);


private:

    // Private Member Functions

        //- Create new streams as required
        void createProbeFiles(const wordList& fieldNames);

        //- Write field values
        template<class Type>
        void writeValues
        (
            const word& fieldName,
            const Field<Type>& values,
            const scalar timeValue
        );

        //- Sample and store/write all applicable sampled fields
        template<class GeoField>
        void performAction
        (
            const fieldGroup<GeoField>& fieldNames,  /* must be sorted */
            unsigned request
        );

        //- Perform sampling action with store/write
        bool performAction(unsigned request);


public:

    // Constructors

        //- Construct from Time and dictionary
        ProbesBase
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );


    //- Destructor
    virtual ~ProbesBase() = default;


    // Member Functions

        //- Enable/disable verbose output
        //  \return old value
        bool verbose(const bool on) noexcept;

        //- Return names of fields to probe
        virtual const wordRes& fieldNames() const noexcept
        {
            return fieldSelection_;
        }

        //- Return const reference to the point prober
        const Prober& prober() const noexcept { return prober_; }

        //- Read the settings from the dictionary
        virtual bool read(const dictionary&);

        //- Sample and store result if the sampleOnExecute is enabled.
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for changes of mesh due to readUpdate
        virtual void readUpdate(const polyMesh::readUpdateState state)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ProbesBase.C"
    #include "ProbesBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
