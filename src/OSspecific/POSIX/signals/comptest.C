// The file sigFpe.C contains code that will not compile if
// - The OS is Linux -> __linux__ is defined
// - The compiler is gcc -> __GNUX__ is defined
// - The libc in use is not glibc, but e.g. musl.
// You can create such a platform by taking a clean alpine linux and installing these packages:
// apk add git meson bash g++ zlib-dev fftw-dev openmpi-dev boost-dev flex-dev cgal-dev
// On this alpine machine, sigFpe.C will not compile since feenableexcept is not defined.
// There is no way to make it compile without touching the source code of sigFpe.C. I'm not going to fix sigFpe (at least not right now). Now if you cannot build OpenFOAM on a machine, it is better to fail a few seconds after the call to `meson setup` than to fail hours after the call to `ninja`. So, I extracted the critical part of sigFpe.C into this file, comptest.C. openfoam/meson.build will check if this file compiles and refuse to build openfoam is this file cannot be compiled.

#if defined(__linux__) && defined(__GNUC__)
    #ifndef __USE_GNU
        #define __USE_GNU      // To use feenableexcept()
    #endif
    #include <fenv.h>
    #include <malloc.h>
#endif

#ifdef __APPLE__
    #include "feexceptErsatz.H"
#endif


#if (defined(__linux__) && defined(__GNUC__)) || defined(__APPLE__)
auto funcptr = feenableexcept;
#endif
