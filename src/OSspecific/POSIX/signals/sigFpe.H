/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2025 OpenCFD Ltd.
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
    Foam::sigFpe

Description
    Set up trapping for floating point exceptions (signal FPE).

    Defined by controlDict InfoSwitch entries:
    - \par trapFpe
      Enable floating point exception trapping.

    - \par setNaN
      Initialization all malloced memory to NaN.
      Combined with \c trapFpe, this causes usage of uninitialized scalars
      to trigger an abort.

    Environment variables:
      - \par FOAM_SIGFPE (true|false)
        overrides \c trapFpe
      - \par FOAM_SETNAN  (true|false)
        overrides \c setNaN

    Note that trapping can be set/removed through the static member functions
    or through the scope of the object (constructor sets trapping; destructor
    restores original). The class behaves as a singleton.

SourceFiles
    sigFpe.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_sigFpe_H
#define Foam_sigFpe_H

#include <algorithm>
#include <limits>
#include <cstddef>  // For std::size_t
#include "scalarFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class UList;

/*---------------------------------------------------------------------------*\
                           Class sigFpe Declaration
\*---------------------------------------------------------------------------*/

class sigFpe
{
    // Private Data

        //- Flag that floating point trapping should be used.
        //  Can override with FOAM_SIGFPE env variable
        static bool switchFpe_;

        //- Flag that NaN initialisation should be used.
        //  Can override with FOAM_SETNAN env variable
        static bool switchNan_;

        //- Floating point trapping currently active?
        static bool sigActive_;

        //- Is NaN memory initialisation currently active?
        static bool nanActive_;


    // Private Member Functions

        //- Handler for caught signals - ends job and prints stack
        static void sigHandler(int);

        //- Fill (float | double) buffer with signalling NaN
        template<class FloatType>
        inline static void fill_with_NaN(FloatType* buf, size_t count)
        {
            if (buf && count)
            {
                const auto val =
                    std::numeric_limits<FloatType>::signaling_NaN();

                // [float|double] so do not need std::uninitialized_fill_n()
                // Can dispatch with
                // - std::execution::par_unseq
                // - std::execution::unseq
                std::fill_n(buf, count, val);
            }
        }


public:

    // Constructors

        //- Constructor calls set() to activate the FPE signal handler if it
        //- was was not previously activate and requested() returns true.
        sigFpe();


    //- Destructor calls unset() to deactivate the FPE signal handler
    //- as required.
    ~sigFpe();


    // Static Member Functions

        //- Check if SIGFPE signals handler is to be enabled.
        //  This is controlled by the trapFpe entry or the FOAM_SIGFPE
        //  environment variable
        static bool requested();

        //- True if SIGFPE handling is currently active.
        static bool active() noexcept { return sigActive_; }

        //- True if NaN memory initialisation is currently active.
        static bool nanActive() noexcept { return nanActive_; }

        //- Set NaN memory initialisation on/off.
        //  \return the previous value
        static bool nanActive(bool on) noexcept
        {
            bool old(nanActive_);
            nanActive_ = on;
            return old;
        }

        //- Activate SIGFPE handler when FOAM_SIGFPE is enabled.
        //- Activate fill memory with signaling_NaN when FOAM_SETNAN is enabled
        static void set(bool verbose=false);

        //- Deactivate SIGFPE handler and NaN memory initialisation
        static void unset(bool verbose=false);

        //- Fill data block with signaling_NaN values if nanActive().
        //- Does a reinterpret to \c Foam::scalar
        static void fillNan_if(void* buf, size_t count)
        {
            if (nanActive_)
            {
                fill_with_NaN
                (
                    reinterpret_cast<Foam::scalar*>(buf),
                    (count/sizeof(Foam::scalar))
                );
            }
        }

        //- Fill data block with signaling_NaN values.
        //- Does a reinterpret to \c Foam::scalar
        static void fillNan(char* buf, size_t count)
        {
            fill_with_NaN
            (
                reinterpret_cast<Foam::scalar*>(buf),
                (count/sizeof(Foam::scalar))
            );
        }

        //- Fill data block with (float) signaling_NaN values
        static void fillNan(UList<float>& list);

        //- Fill data block with (double) signaling_NaN values
        static void fillNan(UList<double>& list);


    // Helpers

        //- Helper to locally ignore SIGFPE handling.
        //  Restores the original state of the SIGFPE handler on destruction.
        class ignore
        {
            //- The signal handler state when entering
            bool wasActive_;

        public:

            //- No copy construct
            ignore(const ignore&) = delete;

            //- No move construct
            ignore(ignore&&) = delete;

            //- No copy assignment
            void operator=(const ignore&) = delete;

            //- No move assignment
            void operator=(ignore&&) = delete;

            //- Constructor deactivates any previously active SIGFPE handler
            ignore()
            :
                wasActive_(Foam::sigFpe::active())
            {
                if (wasActive_)
                {
                    Foam::sigFpe::unset();
                }
            }

            //- Destructor restores the original state of SIGFPE handler
            ~ignore() { reset(); }

            //- Restore the original state of SIGFPE handler
            void reset()
            {
                if (wasActive_)
                {
                    wasActive_ = false;
                    Foam::sigFpe::set();
                }
            }

            //- Same as reset()
            void restore() { reset(); }
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
