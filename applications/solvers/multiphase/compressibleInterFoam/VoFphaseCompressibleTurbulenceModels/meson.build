srcfiles = [
    lnInclude_hack,
    'VoFphaseCompressibleTurbulenceModels.C',
    'compressibleInterPhaseTransportModel.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_twoPhaseMixtureThermo,
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_specie,
    lib_twoPhaseMixture,
    lib_twoPhaseProperties,
    lib_interfaceProperties,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
    lib_incompressibleTransportModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/compressibleInterFoam/VoFphaseCompressibleTurbulenceModels',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/compressibleInterFoam/twoPhaseMixtureThermo',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/phaseCompressible',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/compressibleInterFoam/VoFphaseCompressibleTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_VoFphaseCompressibleTurbulenceModels = library(
    'VoFphaseCompressibleTurbulenceModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_VoFphaseCompressibleTurbulenceModels)
