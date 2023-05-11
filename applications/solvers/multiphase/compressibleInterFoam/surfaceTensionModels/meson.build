srcfiles = [
    lnInclude_hack,
    'liquidProperties/liquidPropertiesSurfaceTension.C',
]
link_with = [
    lib_interfaceProperties,
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_specie,
    lib_thermophysicalProperties,
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/compressibleInterFoam/surfaceTensionModels',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/thermophysicalProperties',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/compressibleInterFoam/surfaceTensionModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_twoPhaseSurfaceTension = library(
    'twoPhaseSurfaceTension',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_twoPhaseSurfaceTension)
