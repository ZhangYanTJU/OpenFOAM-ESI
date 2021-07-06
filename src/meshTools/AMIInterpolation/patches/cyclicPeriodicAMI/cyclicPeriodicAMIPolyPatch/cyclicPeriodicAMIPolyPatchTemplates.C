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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "transformField.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicPeriodicAMIPolyPatch::interpolate
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    // Get the periodic patch
    const coupledPolyPatch& pp
    (
        refCast<const coupledPolyPatch>
        (
            boundaryMesh()[periodicPatchID()]
        )
    );

    if (debug & 2)
    {
        Pout<< "For patch " << this->name()
            << " found periodic:" << pp.name()
            //<< " own:" << pp.owner()
            << " rank:" << pTraits<Type>::rank
            << " parallel:" << pp.parallel()
            //<< " forwardT:" << pp.forwardT()
            //<< " reverseT:" << pp.reverseT()
            << endl;
    }

    if (pTraits<Type>::rank == 0 || pp.parallel())
    {
        // No need for rotation. Pass up to default AMI interpolation.
        return cyclicAMIPolyPatch::interpolate(fld, defaultValues);
    }
    else
    {
        // Transform
        // ~~~~~~~~~
        // - transform to cylindrical coords
        // - do AMI interpolation
        // - transform back to cartesian

        const cyclicAMIPolyPatch& nbrPp = this->neighbPatch();

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

        vector axis(Zero);
        point axisPoint(Zero);
        if (isA<cyclicPolyPatch>(pp))
        {
            axis = refCast<const cyclicPolyPatch>(pp).rotationAxis();
            axisPoint = refCast<const cyclicPolyPatch>(pp).rotationCentre();
        }
        else if (isA<cyclicAMIPolyPatch>(pp))
        {
            axis = refCast<const cyclicAMIPolyPatch>(pp).rotationAxis();
            axisPoint = refCast<const cyclicAMIPolyPatch>(pp).rotationCentre();
        }
        else
        {
            FatalErrorInFunction << "On patch " << name()
                << " have unsupported periodicPatch " << pp.name()
                << exit(FatalError);
        }

        const coordSystem::cylindrical cs(axisPoint, axis);

        auto tlocalFld(tmp<Field<Type>>::New(fld.size()));
        auto& localFld = tlocalFld.ref();
        List<Type> localDeflt(defaultValues.size());

        // Transform to cylindrical coords
        {
            tmp<tensorField> nbrT(cs.R(nbrPp.faceCentres()));
            localFld = Foam::invTransform(nbrT, fld);
            if (defaultValues.size())
            {
                // We get in UList (why? Copied from cyclicAMI). Convert to
                // Field so we can use transformField routines.
                const SubField<Type> defaultSubFld(defaultValues);
                const Field<Type>& defaultFld(defaultSubFld);
                localDeflt = Foam::invTransform(nbrT, defaultFld);
            }
        }


        if (debug&2)
        {
            const pointField& nbrFc = nbrPp.faceCentres();

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

            if (defaultValues.size())
            {
                forAll(defaultValues, i)
                {
                    Pout<< "Defaults At:" << nbrFc[i] << nl
                        << "    cart:" << defaultValues[i] << nl
                        << "    cyli:" << localDeflt[i] << nl
                        << endl;
                }
            }
        }


        // Do the actual interpolation
        auto tinterp = cyclicAMIPolyPatch::interpolate(localFld, localDeflt);

        // Transform back. Result is now at *this
        List<Type> oldResult;
        if (debug&2)
        {
            oldResult = tinterp();
        }

        const pointField& fc = this->faceCentres();
        tmp<tensorField> T(cs.R(fc));
        auto tresult = Foam::transform(T, tinterp);

        if (debug&2)
        {
            const pointField& fc = this->faceCentres();
            forAll(fc, i)
            {
                Pout<< "Results: At local :" << fc[i] << nl
                    << "    cyli:" << oldResult[i] << nl
                    << "    cart:" << tresult()[i] << nl
                    << endl;
            }
        }

        return tresult;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::cyclicPeriodicAMIPolyPatch::interpolate
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    // Get the periodic patch
    const coupledPolyPatch& pp
    (
        refCast<const coupledPolyPatch>
        (
            boundaryMesh()[periodicPatchID()]
        )
    );

    if (pTraits<Type>::rank == 0 || pp.parallel())
    {
        // No need for rotation. Pass up to default AMI interpolation.
        return cyclicAMIPolyPatch::interpolate(tFld, defaultValues);
    }
    else
    {
        return interpolate(tFld(), defaultValues);
    }
}


// ************************************************************************* //
