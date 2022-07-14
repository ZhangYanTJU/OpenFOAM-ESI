/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fullIncompressibleVars::storeOldTimes
(
    GeometricField<Type, PatchField, GeoMesh>& field,
    PtrList<compressedGeometricField<Type, PatchField, GeoMesh>>&
        compressedField,
    const label pos
)
{
    // This object is constructed to store the primal solution of time-step i.
    // The implication is that this happens exactly before solving for i+1,
    // meaning that the time has shifted and so have the instanteneous primal
    // fields. So by the time the function is called, UInst() contains the
    // solution of ith time-step, as an initialization to the solution that it
    // will be computed at the (i+1)th time-step, UInst().oldTime() contains
    // also the solution of the ith time-step and UInst().oldTime().oldTime(),
    // if set, the solution of the (i-1)th time-step. The solution of the
    // (i-2)th time-step is not necessary to be stored for the continuation of
    // the run. So, UInst() and UInst().oldTime().oldTime(), if set, are
    // stored.
    label nOldTimes = min(field.nOldTimes(), a_);
    if (nOldTimes > 1)
    {
        GeometricField<Type, PatchField, GeoMesh>* oldField = &field.oldTime();
        compressedField.setSize(nOldTimes-1);
        for (label i=0; i<nOldTimes-1; i++)
        {
            oldField = &oldField->oldTime();
            compressedField.set
            (
                i,
                compressedGeometricField<Type, PatchField, GeoMesh>::New
                (
                    *oldField,
                    storageParams_,
                    pos,
                    k_
                )
            );
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fullIncompressibleVars::decompressAll
(
    GeometricField<Type, PatchField, GeoMesh>& GeomField,
    compressedGeometricField<Type, PatchField, GeoMesh>& current,
    PtrList<compressedGeometricField<Type, PatchField, GeoMesh>>& oldFields
)
{
    // Due to the implication described in storeOldTimes(...) function, the
    // reference to the GeometricField of the solver that each element of the
    // oldFields_ is not correct. For example, for a second order time-scheme,
    // oldFields_ has only one elememt which holds a reference to the
    // oldTime().oldTime() of the GeometricField known to the solver.
    // Nevertheless, its stored values must be appended to the oldTime() of the
    // GeometricField known to the solver.
    current.decompress();
    if ( GeomField.nOldTimes() > 1 )
    {
        GeometricField<Type, PatchField, GeoMesh>* oldGeomField = &GeomField;
        //- Scanning only what it has been stored:
        //   i) At the first time-step, nOldTimes() may be zero, but upon
        //      reconstruction the noldTimes of the variablesSet of the solver
        //      will be 1 or 2.
        //  ii) if (a=1) the oldTimes are not stored and as a result can not
        //      be retrieved.
        forAll(oldFields, iPtr)
        {
            oldGeomField = &oldGeomField->oldTime();
            oldFields[iPtr].decompress(*oldGeomField);
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fullIncompressibleVars::addStorageMetricsContribution
(
    const compressedGeometricField<Type, PatchField, GeoMesh>& compField,
    const PtrList<compressedGeometricField<Type, PatchField, GeoMesh>>&
        compOldField
)
{
    forAll (storageMetrics_, Ic)
    {
        storageMetrics_[Ic] += compField.storageMetrics()[Ic];
        forAll(compOldField, iPtr)
        {
            storageMetrics_[Ic] += compOldField[iPtr].storageMetrics()[Ic];
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fullIncompressibleVars::calculateAndWrite
(
    const compressedGeometricField<Type, PatchField, GeoMesh>& compField,
    const PtrList<compressedGeometricField<Type, PatchField, GeoMesh>>&
        compOldField,
    const word& name,
    label& i
)
{
    addStorageMetricsContribution(compField, compOldField);
    write(++i);
    writeLog(name);
    storageMetrics_ = Zero;
}


// ************************************************************************* //
