srcfiles = [
    lnInclude_hack,
    'implicitFunctions/implicitFunction.C',
    'implicitFunctions/sphere/sphereImplicitFunction.C',
    'implicitFunctions/sin/sinImplicitFunction.C',
    'implicitFunctions/ellipsoid/ellipsoidImplicitFunction.C',
    'implicitFunctions/paraboloid/paraboloidImplicitFunction.C',
    'implicitFunctions/plane/planeImplicitFunction.C',
    'implicitFunctions/cylinder/cylinderImplicitFunction.C',
    'implicitFunctions/composedFunction/composedFunctionImplicitFunction.C',
]
link_with = []
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/preProcessing/setAlphaField/alphaFieldFunctions',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/setAlphaField/alphaFieldFunctions',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_alphaFieldFunctions = library(
    'alphaFieldFunctions',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_alphaFieldFunctions)
