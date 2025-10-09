/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2025 OpenCFD Ltd.
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
    Foam::regionModels::regionFaModel

Description
    Base class for area region models.

Usage
    Example of the model specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        region          <word>;
        active          <bool>;

        <model>Coeffs
        {
            // subdictionary entries
        }

        // Optional entries
        infoOutput      <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type  | Reqd | Deflt
      region     | Name of operand region              | word  | yes  | -
      area       | Name of the finite-area mesh        | word  | no   | region0
      active     | Flag to activate the model          | bool  | yes  | -
      suffixing  | Suffix hint for model variables     | word  | no   | -
      infoOutput | Flag to activate information output | bool  | no   | false
    \endtable

Note
    The \c suffixing parameter is a hint that individual models may use
    when automatically generating variable names internally. For example,
    a model \em may provide ("T" + suffixing) and ("q" + suffixing)
    as its default names temperature and flux fields, respectively.
    If the user specifies \c suffixing = "_foo", those \em default names
    would then become "T_foo" and "q_foo", respectively.
    Suffixing (true|false|none|default|...) currently selects between
    no suffix and ('_'+region).

SourceFiles
    regionFaModelI.H
    regionFaModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_regionFaModel_H
#define Foam_regionFaModel_H

#include "volMesh.H"
#include "volSurfaceMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                         Class regionFaModel Declaration
\*---------------------------------------------------------------------------*/

class regionFaModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Construct region mesh and fields
        void constructMeshObjects(const dictionary&);

        //- Initialise the region
        void initialise();

        //- Read control parameters from dictionary
        bool init(const dictionary& dict);


protected:

    // Protected Data

        //- Reference to the primary mesh database
        const fvMesh& primaryMesh_;

        //- Reference to the time database
        const Time& time_;

        //- Active flag
        Switch active_;

        //- Active information output
        Switch infoOutput_;

        //- Model name
        const word modelName_;

        //- Suffix hint for automatic model variable names (default: "")
        word suffixHint_;

        //- The finite-area mesh name (default: region0)
        word areaName_;

        //- Region name
        word regionName_;

        //- Model coefficients dictionary
        dictionary coeffs_;

        //- Pointer to the region mesh database
        autoPtr<faMesh> regionMeshPtr_;

        //- Dictionary of output properties
        autoPtr<IOdictionary> outputPropertiesPtr_;

        //- Volume/surface mapping
        autoPtr<volSurfaceMapping> vsmPtr_;


public:

    //- Runtime type information
    TypeName("regionFaModel");


    //- Default name regionFaModel
    static const word regionFaModelName;

    // Constructors

        //- Construct from mesh and name and dict
        regionFaModel
        (
            const fvMesh& mesh,
            const word& regionType,
            const word& modelName,
            const dictionary& dict,
            bool readFields = true
        );

        //- No copy construct
        regionFaModel(const regionFaModel&) = delete;

        //- No copy assignment
        void operator=(const regionFaModel&) = delete;


    //- Destructor
    virtual ~regionFaModel() = default;


    // Member Functions

    // Access

        //- Return the reference to the primary mesh database
        const fvMesh& primaryMesh() const noexcept { return primaryMesh_; }

        //- Return the reference to the time database
        const Time& time() const noexcept { return time_; }

        //- Return the active flag
        Switch active() const noexcept { return active_; }

        //- Return the information flag
        Switch infoOutput() const noexcept { return infoOutput_; }

        //- The finite-area mesh name (extracted from dictionary)
        const word& areaName() const noexcept { return areaName_; }

        //- Return the model name
        const word& modelName() const noexcept { return modelName_; }

        //- Return the region mesh database
        inline const faMesh& regionMesh() const;

        //- Return the region mesh database for manipulation
        inline faMesh& regionMesh();

        //- Return the model coefficients dictionary
        const dictionary& coeffs() const noexcept { return coeffs_; }

        //- Return const access to the output properties dictionary
        inline const IOdictionary& outputProperties() const;

        //- Return output properties dictionary
        inline IOdictionary& outputProperties();

        //- Return the solution dictionary
        inline const dictionary& solution() const;


    // Addressing

        //- List of patch IDs on the primary region coupled to this region
        inline const labelList& primaryPatchIDs() const;

        //- True if patchi on the primary region is coupled to this region
        inline bool isRegionPatch(const label patchi) const;


    // Helper Functions

        //- Return mapping between surface and volume fields
        const volSurfaceMapping& vsm() const;

        //- The suffix hint for automatic model variable names
        const word& suffixHint() const noexcept { return suffixHint_; }

        //- Return the concatenation of \p base and the suffix hint
        inline word suffixed(const char* base) const;

        //- Return the concatenation of \p base and the suffix hint
        inline word suffixed(const std::string& base) const;


    // Evolution

        //- Main driver routing to evolve the region - calls other evolves
        virtual void evolve();

        //- Pre-evolve region
        virtual void preEvolveRegion();

        //- Evolve the region
        virtual void evolveRegion();

        //- Post-evolve region
        virtual void postEvolveRegion();

        //- Courant number of the region
        virtual scalar CourantNumber() const;


    // IO

        //- Provide some feedback
        virtual void info() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionFaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionFaModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
