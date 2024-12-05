/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2025 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolateUntransformed
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    if (owner())
    {
        return AMI().interpolateToSource(fld, defaultValues);
    }
    else
    {
        return neighbPatch().AMI().interpolateToTarget(fld, defaultValues);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolate
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    // Can rotate fields (vector and non-spherical tensors)
    constexpr bool transform_supported = is_rotational_vectorspace_v<Type>;

    [[maybe_unused]]
    autoPtr<coordSystem::cylindrical> cs;

    // Similar to doTransform.

    if constexpr (transform_supported)
    {
        cs.reset(cylindricalCS());
    }

    if (!cs)
    {
        return interpolateUntransformed(fld, defaultValues);
    }

    if constexpr (transform_supported)
    {
        const cyclicAMIPolyPatch& nbrPp = this->neighbPatch();

        if (debug)
        {
            Pout<< "cyclicAMIPolyPatch::interpolate :"
                << " patch:" << this->name()
                << " size:" << this->size()
                << " nbrPatch:" << nbrPp.name()
                << " size:" << nbrPp.size()
                << endl;
        }

        if (fld.size() != nbrPp.size())
        {
            FatalErrorInFunction
                << "Patch:" << this->name()
                << " size:" << this->size()
                << " neighbour patch:" << nbrPp.name()
                << " size:" << nbrPp.size()
                << " fld size:" << fld.size()
                << exit(FatalError);
        }


        Field<Type> localFld(fld.size());

        // Transform to cylindrical coords
        {
            const tensorField nbrT(cs().R(nbrPp.faceCentres()));
            Foam::invTransform(localFld, nbrT, fld);
        }

        if (debug&2)
        {
            const vectorField::subField nbrFc(nbrPp.faceCentres());

            Pout<< "On patch:" << this->name()
                << " size:" << this->size()
                << " fc:" << gAverage(this->faceCentres())
                << " getting remote data from:" << nbrPp.name()
                << " size:" << nbrPp.size()
                << " fc:" << gAverage(nbrFc)
                << endl;

            forAll(fld, i)
            {
                Pout<< "At:" << nbrFc[i] << nl
                    << "    cart:" << fld[i] << nl
                    << "    cyli:" << localFld[i] << nl
                    << endl;
            }
        }


        const tensorField ownT(cs().R(this->faceCentres()));

        Field<Type> localDeflt(defaultValues.size());
        if (defaultValues.size() == size())
        {
            // Transform default values into cylindrical coords (using
            // *this faceCentres)
            // We get in UList (why? Copied from cyclicAMI). Convert to
            // Field so we can use transformField routines.
            const SubField<Type> defaultSubFld(defaultValues);
            const Field<Type>& defaultFld(defaultSubFld);
            Foam::invTransform(localDeflt, ownT, defaultFld);
        }

        // Do the actual interpolation and interpolate back to cartesian
        return Foam::transform
        (
            ownT,
            interpolateUntransformed(localFld, localDeflt)
        );
    }
    else  // (!transform_supported)
    {
        FatalErrorInFunction
            << "CODING ERROR??" << nl
            << "calculated cylindrical coordinate system,"
               " but does not appear to be a vector-space type" << endl
            << Foam::abort(FatalError);
        return nullptr;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolate
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolate(tFld(), defaultValues);
}


template<class Type>
void Foam::cyclicAMIPolyPatch::initInterpolateUntransformed
(
    const Field<Type>& fld,
    labelRange& sendRequests,
    PtrList<List<Type>>& sendBuffers,
    labelRange& recvRequests,
    PtrList<List<Type>>& recvBuffers
) const
{
    const auto& AMI = (owner() ? this->AMI() : neighbPatch().AMI());

    if (AMI.distributed() && AMI.comm() != -1)
    {
        const auto& map = (owner() ? AMI.tgtMap() : AMI.srcMap());

        // Insert send/receive requests (non-blocking)
        map.send
        (
            fld,
            sendRequests,
            sendBuffers,
            recvRequests,
            recvBuffers,
            3894+this->index()  // unique offset + patch index
        );
    }
}


template<class Type>
void Foam::cyclicAMIPolyPatch::initInterpolate
(
    const Field<Type>& fld,
    labelRange& sendRequests,
    PtrList<List<Type>>& sendBuffers,
    labelRange& recvRequests,
    PtrList<List<Type>>& recvBuffers
) const
{
    const auto& AMI = (owner() ? this->AMI() : neighbPatch().AMI());

    if (!AMI.distributed() || AMI.comm() == -1)
    {
        return;
    }

    // Can rotate fields (vector and non-spherical tensors)
    constexpr bool transform_supported = is_rotational_vectorspace_v<Type>;

    [[maybe_unused]]
    autoPtr<coordSystem::cylindrical> cs;

    if constexpr (transform_supported)
    {
        cs.reset(cylindricalCS());
    }

    if (!cs)
    {
        initInterpolateUntransformed
        (
            fld,
            sendRequests,
            sendBuffers,
            recvRequests,
            recvBuffers
        );
    }
    else if constexpr (transform_supported)
    {
        const cyclicAMIPolyPatch& nbrPp = this->neighbPatch();

        Field<Type> localFld(fld.size());

        // Transform to cylindrical coords
        {
            const tensorField nbrT(cs().R(nbrPp.faceCentres()));
            Foam::invTransform(localFld, nbrT, fld);
        }

        initInterpolateUntransformed
        (
            localFld,
            sendRequests,
            sendBuffers,
            recvRequests,
            recvBuffers
        );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolate
(
    const Field<Type>& localFld,
    const labelRange& requests,
    const PtrList<List<Type>>& recvBuffers,
    const UList<Type>& defaultValues
) const
{
    auto tresult = tmp<Field<Type>>::New(this->size(), Zero);

    const auto& AMI = (owner() ? this->AMI() : neighbPatch().AMI());

    Field<Type> work;
    if (AMI.distributed())
    {
        if (AMI.comm() == -1)
        {
            return tresult;
        }

        const auto& map = (owner() ? AMI.tgtMap() : AMI.srcMap());

        // Receive (= copy) data from buffers into work. TBD: receive directly
        // into slices of work.
        map.receive
        (
            requests,
            recvBuffers,
            work,
            3894+this->index()  // unique offset + patch index
        );
    }
    const Field<Type>& fld = (AMI.distributed() ? work : localFld);

    // Rotate fields (vector and non-spherical tensors)
    constexpr bool transform_supported = is_rotational_vectorspace_v<Type>;

    [[maybe_unused]]
    autoPtr<coordSystem::cylindrical> cs;

    // Rotate fields (vector and non-spherical tensors)
    if constexpr (transform_supported)
    {
        cs.reset(cylindricalCS());
    }

    if (!cs)
    {
        AMI.weightedSum
        (
            owner(),
            fld,
            tresult.ref(),
            defaultValues
        );
    }
    else if constexpr (transform_supported)
    {
        const tensorField ownT(cs().R(this->faceCentres()));

        Field<Type> localDeflt(defaultValues.size());
        if (defaultValues.size() == size())
        {
            // Transform default values into cylindrical coords (using
            // *this faceCentres)
            // We get in UList (why? Copied from cyclicAMI). Convert to
            // Field so we can use transformField routines.
            const SubField<Type> defaultSubFld(defaultValues);
            const Field<Type>& defaultFld(defaultSubFld);
            Foam::invTransform(localDeflt, ownT, defaultFld);
        }

        AMI.weightedSum
        (
            owner(),
            fld,
            tresult.ref(),
            localDeflt
        );

        // Transform back
        Foam::transform(tresult.ref(), ownT, tresult());
    }

    return tresult;
}


template<class Type, class CombineOp>
void Foam::cyclicAMIPolyPatch::interpolate
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    //- Commented out for now since called with non-primitives (e.g. wallPoint
    //  from FaceCellWave) - missing Foam::transform, Foam::invTransform
    /*
     *
    // Rotate fields (vector and non-spherical tensors)
    constexpr bool transform_supported = is_rotational_vectorspace_v<Type>;

    [[maybe_unused]]
    autoPtr<coordSystem::cylindrical> cs;

    // Rotate fields (vector and non-spherical tensors)
    if constexpr (transform_supported)
    {
        cs.reset(cylindricalCS());
    }

    if (cs)
    {
        const cyclicAMIPolyPatch& nbrPp = this->neighbPatch();

        // Transform to cylindrical coords
        {
            const tensorField nbrT(cs().R(nbrPp.faceCentres()));
            Foam::invTransform(result, nbrT, result);
        }

        const tensorField ownT(cs().R(this->faceCentres()));

        Field<Type> localDeflt(defaultValues.size());
        if (defaultValues.size() == size())
        {
            // Transform default values into cylindrical coords (using
            // *this faceCentres)
            // We get in UList (why? Copied from cyclicAMI). Convert to
            // Field so we can use transformField routines.
            const SubField<Type> defaultSubFld(defaultValues);
            const Field<Type>& defaultFld(defaultSubFld);
            Foam::invTransform(localDeflt, ownT, defaultFld);
        }

        // Do actual AMI interpolation
        if (owner())
        {
            AMI().interpolateToSource
            (
                fld,
                cop,
                result,
                localDeflt
            );
        }
        else
        {
            neighbPatch().AMI().interpolateToTarget
            (
                fld,
                cop,
                result,
                localDeflt
            );
        }

        // Transform back. Result is now at *this
        Foam::transform(result, ownT, result);
    }
    else
    */
    {
        if (owner())
        {
            AMI().interpolateToSource
            (
                fld,
                cop,
                result,
                defaultValues
            );
        }
        else
        {
            neighbPatch().AMI().interpolateToTarget
            (
                fld,
                cop,
                result,
                defaultValues
            );
        }
    }
}


// ************************************************************************* //
