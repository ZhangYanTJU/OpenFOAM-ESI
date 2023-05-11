srcfiles = [
    lnInclude_hack,
    'DPMTurbulenceModels.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_compressibleTransportModels,
    lib_incompressibleTransportModels,
    lib_turbulenceModels,
    lib_incompressibleTurbulenceModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/lagrangian/DPMFoam/DPMTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + meson.source_root() / 'src/transportModels/incompressible/singlePhaseTransportModel',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/phaseIncompressible',
    '-I' + recursive_include_dirs / 'applications/solvers/lagrangian/DPMFoam/DPMTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_DPMTurbulenceModels = library(
    'DPMTurbulenceModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_DPMTurbulenceModels)
