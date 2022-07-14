/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 PCOpt/NTUA
    Copyright (C) 2022      FOSS GP
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
    Foam::fullIncompressibleVars

Description
    Base class storing the current and the oldTime primal variables instances

\*---------------------------------------------------------------------------*/

#ifndef fullIncompressibleVars_H
#define fullIncompressibleVars_H

#include "compressedIncompressibleVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fullIncompressibleVars Declaration
\*---------------------------------------------------------------------------*/

class fullIncompressibleVars
:
    public compressedIncompressibleVars
{
protected:

    // Protected data

        //- Total number of time solutions to be stored for each time instance.
        //  E.g. for (a=1) only the solution at the current time-step will be
        //  stored, for (a=2) the solution at the current and the previous
        //  tim-step will be stored.
        label a_;

        //- PtrLists containing the oldTimes of all primal fields.
        //  The fields of the current time-step are stored upon construction of
        //  the compressedGeometricField class.
        PtrList<compressedVolScalarField> pOld_;
        PtrList<compressedSurfaceScalarField> phiOld_;
        PtrList<compressedVolVectorField> UOld_;
        PtrList<compressedVolScalarField> TMVar1Old_;
        PtrList<compressedVolScalarField> TMVar2Old_;


    // Protected Member Functions

        //- Store oldTimes (if any) of a
        //- GeometricField<Type, PatchField, GeoMesh>
        template<class Type, template<class> class PatchField, class GeoMesh>
        void storeOldTimes
        (
            GeometricField<Type, PatchField, GeoMesh>& field,
            PtrList<compressedGeometricField<Type, PatchField, GeoMesh>>&
                compressedField,
            const label pos
        );

        //- Decompress field and its old times into a GeometricField
        template<class Type, template<class> class PatchField, class GeoMesh>
        void decompressAll
        (
            GeometricField<Type, PatchField, GeoMesh>& field,
            compressedGeometricField<Type, PatchField, GeoMesh>& current,
            PtrList<compressedGeometricField<Type, PatchField, GeoMesh>>&
                oldFields
        );

        //- Calculate and write compression metrics for a
        //- compressedGeometricField
        template<class Type, template<class> class PatchField, class GeoMesh>
        void addStorageMetricsContribution
        (
            const compressedGeometricField<Type, PatchField, GeoMesh>&
                compField,
            const PtrList<compressedGeometricField<Type, PatchField, GeoMesh>>&
                compOldField
        );

        //- Calculate and write compression metrics for a
        //- compressedGeometricField
        template<class Type, template<class> class PatchField, class GeoMesh>
        void calculateAndWrite
        (
            const compressedGeometricField<Type, PatchField, GeoMesh>&
                compField,
            const PtrList<compressedGeometricField<Type, PatchField, GeoMesh>>&
                compOldField,
            const word& name,
            label& i
        );

        //- Function to calculate and write the storage metrics of the
        //- checkpoint
        virtual void calculateStorageMetrics() override;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fullIncompressibleVars(const fullIncompressibleVars&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const fullIncompressibleVars&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("fullIncompressibleVars");


    // Constructors

        fullIncompressibleVars
        (
            incompressibleVars& vs,
            storageParameters& storageParams,
            label a = 2
        );


    //- Destructor
    virtual ~fullIncompressibleVars() = default;


    // Member Functions

        //- Function to compress the incompressibleVars
        void compress() override;

        //- Function to decompress the incompressible vars and
        //  update the incompressibleVars with the decompressed fields
        void decompress
        (
            incompressibleVars& vars
        ) override;

        void decompress();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fullIncompressibleVarsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
