/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2025 OpenCFD Ltd.
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

#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "UPstreamWrapping.H"

#include <cinttypes>

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Special reductions for bool

void Foam::UPstream::reduceAnd(bool& value, const label comm)
{
    PstreamDetail::allReduce(&value, 1, MPI_C_BOOL, MPI_LAND, comm);
}


void Foam::UPstream::reduceOr(bool& value, const label comm)
{
    PstreamDetail::allReduce(&value, 1, MPI_C_BOOL, MPI_LOR, comm);
}


void Foam::reduce
(
    bool& value,
    const andOp<bool>&,
    const int tag,  /* (unused) */
    const label comm
)
{
    PstreamDetail::allReduce(&value, 1, MPI_C_BOOL, MPI_LAND, comm);
}


void Foam::reduce
(
    bool& value,
    const orOp<bool>&,
    const int tag,  /* (unused) */
    const label comm
)
{
    PstreamDetail::allReduce(&value, 1, MPI_C_BOOL, MPI_LOR, comm);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static inline bool is_basic_dataType(Foam::UPstream::dataTypes id) noexcept
{
    return
    (
        int(id) >= int(Foam::UPstream::dataTypes::Basic_begin)
     && int(id)  < int(Foam::UPstream::dataTypes::Basic_end)
    );
}

static inline bool is_reduce_opCode(Foam::UPstream::opCodes id) noexcept
{
    return
    (
        int(id) >= int(Foam::UPstream::opCodes::Basic_begin)
     && int(id)  < int(Foam::UPstream::opCodes::Basic_end)
    );
}


namespace
{

using namespace Foam;

// Local function to print some error information
void printErrorMessage
(
    const void* values,
    const UPstream::dataTypes datatype_id,
    const UPstream::opCodes opcode_id
)
{
    FatalError
        << "Bad input for reduce(): likely a programming problem\n";

    if (!is_basic_dataType(datatype_id))
    {
        FatalError<< "    Non-basic data tyoe (" << int(datatype_id) << ")\n";
    }

    if (!is_reduce_opCode(opcode_id))
    {
        FatalError<< "    Invalid reduce op (" << int(opcode_id) << ")\n";
    }

    if (values == nullptr)
    {
        FatalError<< "    nullptr for values\n";
    }
    FatalError<< Foam::endl;
}

} // End anonymous namespace


// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

void Foam::UPstream::mpi_reduce
(
    void* values,                          // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const UPstream::opCodes opCodeId,      // Proper code passed by caller
    const int communicator,                // Index into MPICommunicators_
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
    MPI_Op optype = PstreamGlobals::getOpCode(opCodeId);

    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do - ignore
        return;
    }
    if
    (
        FOAM_UNLIKELY
        (
            !is_basic_dataType(dataTypeId)
         || !is_reduce_opCode(opCodeId)
         || (values == nullptr)
        )
    )
    {
        FatalErrorInFunction;
        printErrorMessage(values, dataTypeId, opCodeId);
        FatalError << Foam::abort(FatalError);
    }

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_reduce] : "
            << " op:" << int(opCodeId)
            << " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << Foam::endl;
    }

    {
        // Regular reduce

        PstreamDetail::reduce0
        (
            values,
            count,
            datatype,
            optype,
            communicator,
            req
        );
    }
}


void Foam::UPstream::mpi_allreduce
(
    void* values,                          // Type checking done by caller
    int count,
    const UPstream::dataTypes dataTypeId,  // Proper type passed by caller
    const UPstream::opCodes opCodeId,      // Proper code passed by caller
    const int communicator,                // Index into MPICommunicators_
    UPstream::Request* req
)
{
    MPI_Datatype datatype = PstreamGlobals::getDataType(dataTypeId);
    MPI_Op optype = PstreamGlobals::getOpCode(opCodeId);

    if (!count || !UPstream::is_parallel(communicator))
    {
        // Nothing to do - ignore
        return;
    }
    if
    (
        FOAM_UNLIKELY
        (
            !is_basic_dataType(dataTypeId)
         || !is_reduce_opCode(opCodeId)
         || (values == nullptr)
        )
    )
    {
        FatalErrorInFunction;
        printErrorMessage(values, dataTypeId, opCodeId);
        FatalError << Foam::abort(FatalError);
    }

    if (FOAM_UNLIKELY(UPstream::debug))
    {
        Perr<< "[mpi_allreduce] : "
            << " op:" << int(opCodeId)
            << " type:" << int(dataTypeId) << " count:" << count
            << " comm:" << communicator
            << Foam::endl;
    }

    {
        // Regular allReduce

        PstreamDetail::allReduce
        (
            values,
            count,
            datatype,
            optype,
            communicator,
            req
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Common reductions

#undef  Pstream_CommonReductions
#define Pstream_CommonReductions(Native, TaggedType)                          \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const minOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        values, size, TaggedType, MPI_MIN, comm                               \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const maxOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        values, size, TaggedType, MPI_MAX, comm                               \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        values, size, TaggedType, MPI_SUM, comm                               \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const minOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        &value, 1, TaggedType, MPI_MIN, comm                                  \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const maxOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        &value, 1, TaggedType, MPI_MAX, comm                                  \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        &value, 1, TaggedType, MPI_SUM, comm                                  \
    );                                                                        \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Floating-point reductions

#undef  Pstream_FloatReductions
#define Pstream_FloatReductions(Native, TaggedType)                           \
                                                                              \
Pstream_CommonReductions(Native, TaggedType);                                 \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const sumOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm,                                                         \
    UPstream::Request& req                                                    \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        values, size, TaggedType, MPI_SUM, comm, &req                         \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const sumOp<Native>&,                                                     \
    const int tag,  /* (unused) */                                            \
    const label comm,                                                         \
    UPstream::Request& req                                                    \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        &value, 1, TaggedType, MPI_SUM, comm, &req                            \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::sumReduce                                                          \
(                                                                             \
    Native& value,                                                            \
    label& count,                                                             \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    if (UPstream::is_parallel(comm))                                          \
    {                                                                         \
        Native values[2];                                                     \
        values[0] = static_cast<Native>(count);                               \
        values[1] = value;                                                    \
                                                                              \
        PstreamDetail::allReduce<Native>                                      \
        (                                                                     \
            values, 2, TaggedType, MPI_SUM, comm                              \
        );                                                                    \
                                                                              \
        count = static_cast<label>(values[0]);                                \
        value = values[1];                                                    \
    }                                                                         \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Bitwise reductions

#undef  Pstream_BitwiseReductions
#define Pstream_BitwiseReductions(Native, TaggedType)                         \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native values[],                                                          \
    const int size,                                                           \
    const bitOrOp<Native>&,                                                   \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        values, size, TaggedType, MPI_BOR, comm                               \
    );                                                                        \
}                                                                             \
                                                                              \
void Foam::reduce                                                             \
(                                                                             \
    Native& value,                                                            \
    const bitOrOp<Native>&,                                                   \
    const int tag,  /* (unused) */                                            \
    const label comm                                                          \
)                                                                             \
{                                                                             \
    PstreamDetail::allReduce<Native>                                          \
    (                                                                         \
        &value, 1, TaggedType, MPI_BOR, comm                                  \
    );                                                                        \
}                                                                             \

\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Pstream_CommonReductions(int32_t, MPI_INT32_T);
Pstream_CommonReductions(int64_t, MPI_INT64_T);
Pstream_CommonReductions(uint32_t, MPI_UINT32_T);
Pstream_CommonReductions(uint64_t, MPI_UINT64_T);

Pstream_FloatReductions(float, MPI_FLOAT);
Pstream_FloatReductions(double, MPI_DOUBLE);

Pstream_BitwiseReductions(unsigned char, MPI_UNSIGNED_CHAR);
Pstream_BitwiseReductions(unsigned int, MPI_UNSIGNED);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef Pstream_CommonReductions
#undef Pstream_FloatReductions
#undef Pstream_BitwiseReductions


// ************************************************************************* //
