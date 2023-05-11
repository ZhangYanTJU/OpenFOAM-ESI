srcfiles = [
    lnInclude_hack,
    'compressibleTwoPhaseMixtureTurbulenceModels.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_incompressibleTransportModels,
    lib_immiscibleIncompressibleTwoPhaseMixture,
    lib_twoPhaseProperties,
    lib_interfaceProperties,
    lib_turbulenceModels,
    lib_incompressibleTurbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/MPPICInterFoam/compressibleTwoPhaseMixtureTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/transportModels/incompressible',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/transportModels/immiscibleIncompressibleTwoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/phaseCompressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/phaseIncompressible',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/MPPICInterFoam/compressibleTwoPhaseMixtureTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_compressibleTwoPhaseMixtureTurbulenceModels = library(
    'compressibleTwoPhaseMixtureTurbulenceModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_compressibleTwoPhaseMixtureTurbulenceModels)
