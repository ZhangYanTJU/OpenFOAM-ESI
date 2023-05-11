srcfiles = [
    lnInclude_hack,
    'compressibleMultiPhaseTurbulenceModels.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_fluidThermophysicalModels,
    lib_specie,
    lib_incompressibleTransportModels,
    lib_compressibleTransportModels,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/phaseSystemModels/multiphaseInter/compressibleMultiPhaseTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/multiphaseInter/phasesSystem',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/phaseCompressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/phaseIncompressible',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/multiphaseInter/compressibleMultiPhaseTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_compressibleMultiPhaseTurbulenceModels = library(
    'compressibleMultiPhaseTurbulenceModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_compressibleMultiPhaseTurbulenceModels)
