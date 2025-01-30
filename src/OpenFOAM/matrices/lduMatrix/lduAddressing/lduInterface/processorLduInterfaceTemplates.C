/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "processorLduInterface.H"
#include "IPstream.H"
#include "OPstream.H"

// * * * * * * * * * * * * * * * Member Functions * * *  * * * * * * * * * * //

template<class Type>
void Foam::processorLduInterface::send
(
    const UPstream::commsTypes commsType,
    const UList<Type>& fld
) const
{
    if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << Foam::abort(FatalError);
    }

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::scheduled
    )
    {
        UOPstream::write
        (
            commsType,
            neighbProcNo(),
            fld,
            tag(),
            comm()
        );
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        if (fld.empty())
        {
            // Can skip empty messages
            return;
        }

        const label nBytes = fld.size_bytes();

        resizeBuf(byteSendBuf_, nBytes);
        resizeBuf(byteRecvBuf_, nBytes);

        std::memcpy
        (
            static_cast<void*>(byteSendBuf_.data()), fld.cdata(), nBytes
        );

        UIPstream::read
        (
            commsType,
            neighbProcNo(),
            byteRecvBuf_.data(),
            nBytes,
            tag(),
            comm()
        );

        UOPstream::write
        (
            commsType,
            neighbProcNo(),
            byteSendBuf_.cdata(),
            nBytes,
            tag(),
            comm()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType)
            << exit(FatalError);
    }
}


template<class Type>
void Foam::processorLduInterface::receive
(
    const UPstream::commsTypes commsType,
    UList<Type>& fld
) const
{
    if constexpr (!is_contiguous_v<Type>)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << Foam::abort(FatalError);
    }

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::scheduled
    )
    {
        UIPstream::read
        (
            commsType,
            neighbProcNo(),
            fld,
            tag(),
            comm()
        );
    }
    else if (commsType == UPstream::commsTypes::nonBlocking)
    {
        const label nBytes = fld.size_bytes();

        std::memcpy
        (
            static_cast<void*>(fld.data()), byteRecvBuf_.cdata(), nBytes
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType)
            << exit(FatalError);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::processorLduInterface::receive
(
    const UPstream::commsTypes commsType,
    const label size
) const
{
    auto tfld = tmp<Field<Type>>::New(size);
    receive(commsType, tfld.ref());
    return tfld;
}


template<class Type>
void Foam::processorLduInterface::compressedSend
(
    const UPstream::commsTypes commsType,
    const UList<Type>& f
) const
{
    if constexpr
    (
        std::is_integral_v<Type>
     || (sizeof(float) == sizeof(scalar))
    )
    {
        // No compression if integral or scalar is float
        this->send(commsType, f);
    }
    else if (f.empty() || !UPstream::floatTransfer)
    {
        // No compression
        this->send(commsType, f);
    }
    else  // (!f.empty() && UPstream::floatTransfer)
    {
        static const label nCmpts = (sizeof(Type)/sizeof(scalar));
        const label nm1 = (f.size() - 1)*nCmpts;
        const label nBytes = f.size()*nCmpts*sizeof(float);

        const scalar *sArray = reinterpret_cast<const scalar*>(f.cdata());
        const scalar *slast = &sArray[nm1];
        resizeBuf(byteSendBuf_, nBytes);
        float *fArray = reinterpret_cast<float*>(byteSendBuf_.data());

        for (label i=0; i<nm1; i++)
        {
            fArray[i] = sArray[i] - slast[i%nCmpts];
        }

        reinterpret_cast<Type&>(fArray[nm1]) = f.last();

        if
        (
            commsType == UPstream::commsTypes::buffered
         || commsType == UPstream::commsTypes::scheduled
        )
        {
            UOPstream::write
            (
                commsType,
                neighbProcNo(),
                byteSendBuf_.cdata(),
                nBytes,
                tag(),
                comm()
            );
        }
        else if (commsType == UPstream::commsTypes::nonBlocking)
        {
            resizeBuf(byteRecvBuf_, nBytes);

            UIPstream::read
            (
                commsType,
                neighbProcNo(),
                byteRecvBuf_.data(),
                nBytes,
                tag(),
                comm()
            );

            UOPstream::write
            (
                commsType,
                neighbProcNo(),
                byteSendBuf_.cdata(),
                nBytes,
                tag(),
                comm()
            );
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported communications type " << int(commsType)
                << exit(FatalError);
        }
    }
}


template<class Type>
void Foam::processorLduInterface::compressedReceive
(
    const UPstream::commsTypes commsType,
    UList<Type>& f
) const
{
    if constexpr
    (
        std::is_integral_v<Type>
     || (sizeof(float) == sizeof(scalar))
    )
    {
        // No compression if integral or scalar is float
        this->receive<Type>(commsType, f);
    }
    else if (f.empty() || !UPstream::floatTransfer)
    {
        // Nothing to compress
        this->receive<Type>(commsType, f);
    }
    else  // (!f.empty() && UPstream::floatTransfer)
    {
        static const label nCmpts = (sizeof(Type)/sizeof(scalar));
        const label nm1 = (f.size() - 1)*nCmpts;
        const label nBytes = f.size()*nCmpts*sizeof(float);

        if
        (
            commsType == UPstream::commsTypes::buffered
         || commsType == UPstream::commsTypes::scheduled
        )
        {
            resizeBuf(byteRecvBuf_, nBytes);

            UIPstream::read
            (
                commsType,
                neighbProcNo(),
                byteRecvBuf_.data(),
                nBytes,
                tag(),
                comm()
            );
        }
        else if (commsType != UPstream::commsTypes::nonBlocking)
        {
            FatalErrorInFunction
                << "Unsupported communications type " << int(commsType)
                << exit(FatalError);
        }

        const float *fArray =
            reinterpret_cast<const float*>(byteRecvBuf_.cdata());
        f.last() = reinterpret_cast<const Type&>(fArray[nm1]);
        scalar *sArray = reinterpret_cast<scalar*>(f.data());
        const scalar *slast = &sArray[nm1];

        for (label i=0; i<nm1; i++)
        {
            sArray[i] = fArray[i] + slast[i%nCmpts];
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::processorLduInterface::compressedReceive
(
    const UPstream::commsTypes commsType,
    const label size
) const
{
    auto tfld = tmp<Field<Type>>::New(size);
    compressedReceive(commsType, tfld.ref());
    return tfld;
}


// ************************************************************************* //
