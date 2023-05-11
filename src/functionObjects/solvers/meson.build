srcfiles = [
    lnInclude_hack,
    'scalarTransport/scalarTransport.C',
    'energyTransport/energyTransport.C',
    'electricPotential/electricPotential.C',
]
link_with = [
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_fluidThermophysicalModels,
    lib_incompressibleTransportModels,
    lib_compressibleTransportModels,
    lib_turbulenceModels,
    lib_incompressibleTurbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_specie,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/functionObjects/solvers',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + meson.source_root() / 'src/transportModels/incompressible/singlePhaseTransportModel',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/functionObjects/solvers',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_solverFunctionObjects = library(
    'solverFunctionObjects',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_solverFunctionObjects)
