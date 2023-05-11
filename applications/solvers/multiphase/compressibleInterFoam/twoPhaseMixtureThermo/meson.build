srcfiles = [
    lnInclude_hack,
    'twoPhaseMixtureThermo.C',
]
link_with = [
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_specie,
    lib_twoPhaseMixture,
    lib_interfaceProperties,
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/compressibleInterFoam/twoPhaseMixtureThermo',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/compressibleInterFoam/twoPhaseMixtureThermo',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_twoPhaseMixtureThermo = library(
    'twoPhaseMixtureThermo',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_twoPhaseMixtureThermo)
