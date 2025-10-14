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
        // Only creates the co-ord system if using periodic AMI
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
        if (defaultValues.size() != 0 && defaultValues.size() == size())
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
    labelRange& recvRequests,
    PtrList<List<Type>>& sendBuffers,
    PtrList<List<Type>>& recvBuffers,

    labelRange& sendRequests1,
    labelRange& recvRequests1,
    PtrList<List<Type>>& sendBuffers1,
    PtrList<List<Type>>& recvBuffers1
) const
{
    const auto& AMI = (owner() ? this->AMI() : neighbPatch().AMI());

    if (AMI.distributed() && AMI.comm() != -1)
    {
        const auto& cache = AMI.cache();

        if (cache.index0() == -1 && cache.index1() == -1)
        {
            // No caching

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
        else
        {
            // Caching is active

            cache.setDirection(owner());

            if (cache.index0() != -1)
            {
                const auto& map0 = cache.cTgtMapPtr0()();

                // Insert send/receive requests (non-blocking)
                map0.send
                (
                    fld,
                    sendRequests,
                    sendBuffers,
                    recvRequests,
                    recvBuffers,
                    3894+this->index()  // unique offset + patch index
                );
            }

            if (cache.index1() != -1)
            {
                const auto& map1 = cache.cTgtMapPtr1()();

                // Insert send/receive requests (non-blocking)
                map1.send
                (
                    fld,
                    sendRequests1,
                    sendBuffers1,
                    recvRequests1,
                    recvBuffers1,
                    3895+this->index()  // unique offset + patch index
                );
            }
        }
    }
}


template<class Type>
void Foam::cyclicAMIPolyPatch::initInterpolate
(
    const Field<Type>& fld,
    labelRange& sendRequests,
    labelRange& recvRequests,
    PtrList<List<Type>>& sendBuffers,
    PtrList<List<Type>>& recvBuffers,

    labelRange& sendRequests1,
    labelRange& recvRequests1,
    PtrList<List<Type>>& sendBuffers1,
    PtrList<List<Type>>& recvBuffers1
) const
{
    const auto& AMI = (owner() ? this->AMI() : neighbPatch().AMI());

    if (!AMI.distributed() || AMI.comm() == -1)
    {
        return;
    }

    // Can rotate fields (vector and non-spherical tensors)
    constexpr bool transform_supported = is_rotational_vectorspace_v<Type>;

    if constexpr (transform_supported)
    {
        // Only creates the co-ord system if using periodic AMI
        // - convert to cylindrical coordinate system
        auto cs = cylindricalCS();

        if (cs)
        {
            Field<Type> localFld(fld.size());
            const cyclicAMIPolyPatch& nbrPp = this->neighbPatch();
            const tensorField nbrT(cs().R(nbrPp.faceCentres()));
            Foam::invTransform(localFld, nbrT, fld);

            initInterpolateUntransformed
            (
                localFld,
                sendRequests,
                recvRequests,
                sendBuffers,
                recvBuffers,

                sendRequests1,
                recvRequests1,
                sendBuffers1,
                recvBuffers1
            );

            return;
        }
    }

    initInterpolateUntransformed
    (
        fld,
        sendRequests,
        recvRequests,
        sendBuffers,
        recvBuffers,

        sendRequests1,
        recvRequests1,
        sendBuffers1,
        recvBuffers1
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicAMIPolyPatch::interpolate
(
    const Field<Type>& localFld,
    const labelRange& requests,
    const PtrList<List<Type>>& recvBuffers,
    const labelRange& requests1,
    const PtrList<List<Type>>& recvBuffers1,
    const UList<Type>& defaultValues
) const
{
    // Note: cannot be localFld.size() -> is set to null for distributed AMI
    auto tresult = tmp<Field<Type>>::New(this->size(), Zero);

    const auto& AMI = (owner() ? this->AMI() : neighbPatch().AMI());

    const auto& cache = AMI.cache();
    cache.setDirection(owner());

    Field<Type> work;
    Field<Type> work1;
    if (AMI.distributed())
    {
        if (AMI.comm() == -1)
        {
            return tresult;
        }

        if (cache.index0() == -1 && cache.index1() == -1)
        {
            // No caching
            const auto& map = (owner() ? AMI.tgtMap() : AMI.srcMap());

            // Receive (= copy) data from buffers into work. TBD: receive
            // directly into slices of work.
            map.receive
            (
                requests,
                recvBuffers,
                work,
                3894+this->index()  // unique offset + patch index
            );
        }
        else
        {
            // Using AMI cache

            if (cache.index0() != -1)
            {
                cache.cTgtMapPtr0()().receive
                (
                    requests,
                    recvBuffers,
                    work,
                    3894+this->index()  // unique offset + patch index
                );
            }

            if (cache.index1() != -1)
            {
                cache.cTgtMapPtr1()().receive
                (
                    requests1,
                    recvBuffers1,
                    work1,
                    3895+this->index()  // unique offset + patch index
                );
            }
        }
    }

    const Field<Type>& fld = (AMI.distributed() ? work : localFld);
    const Field<Type>& fld1 = (AMI.distributed() ? work1 : localFld);

    // Rotate fields (vector and non-spherical tensors)
    constexpr bool transform_supported = is_rotational_vectorspace_v<Type>;

    // Rotate fields (vector and non-spherical tensors) for periodic AMI
    tensorField ownTransform;
    Field<Type> localDeflt;

    if constexpr (transform_supported)
    {
        // Only creates the co-ord system if using periodic AMI
        // - convert to cylindrical coordinate system
        auto cs = cylindricalCS();

        if (cs)
        {
            ownTransform = cs().R(this->faceCentres());
            localDeflt = defaultValues;

            if (defaultValues.size() == size())
            {
                // Transform default values into cylindrical coords (using
                // *this faceCentres)
                // We get in UList (why? Copied from cyclicAMI). Convert to
                // Field so we can use transformField routines.
                const SubField<Type> defaultSubFld(defaultValues);
                const Field<Type>& defaultFld(defaultSubFld);
                Foam::invTransform(localDeflt, ownTransform, defaultFld);
            }
        }
    }

    const auto& localDefaultValues =
        localDeflt.size() ? localDeflt : defaultValues;

    if (cache.index0() == -1 && cache.index1() == -1)
    {
        // No caching
        AMI.weightedSum
        (
            owner(),
            fld,
            tresult.ref(),
            localDefaultValues
        );

        // Transform back
        if (ownTransform.size())
        {
            Foam::transform(tresult.ref(), ownTransform, tresult());
        }

        return tresult;
    }
    else
    {
        if (cache.index0() != -1)
        {
            AMIInterpolation::weightedSum
            (
                AMI.lowWeightCorrection(),
                cache.cSrcAddress0(),
                cache.cSrcWeights0(),
                cache.cSrcWeightsSum0(),
                fld,
                multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>()),
                tresult.ref(),
                localDefaultValues
            );

            if (ownTransform.size())
            {
                Foam::transform(tresult.ref(), ownTransform, tresult());
            }

            // Assuming cache weight is zero when index1 is inactive (==-1)
            tresult.ref() *= (1 - cache.weight());
        }

        if (cache.index1() != -1)
        {
            auto tresult1 = tmp<Field<Type>>::New(this->size(), Zero);

            AMIInterpolation::weightedSum
            (
                AMI.lowWeightCorrection(),
                cache.cSrcAddress1(),
                cache.cSrcWeights1(),
                cache.cSrcWeightsSum1(),
                fld1,
                multiplyWeightedOp<Type, plusEqOp<Type>>(plusEqOp<Type>()),
                tresult1.ref(),
                localDefaultValues
            );

            if (ownTransform.size())
            {
                Foam::transform(tresult1.ref(), ownTransform, tresult1());
            }

            tresult1.ref() *= cache.weight();
            tresult.ref() += tresult1();
        }

        return tresult;
    }
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
