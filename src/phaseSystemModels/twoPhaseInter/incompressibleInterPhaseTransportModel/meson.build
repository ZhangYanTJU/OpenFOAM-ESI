srcfiles = [
    lnInclude_hack,
    'incompressibleInterPhaseTransportModels.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_fluidThermophysicalModels,
    lib_specie,
    lib_incompressibleTransportModels,
    lib_turbulenceModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/phaseSystemModels/twoPhaseInter/incompressibleInterPhaseTransportModel',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/twoPhaseInter/VoFphaseIncompressibleTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/immiscibleIncompressibleTwoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/transportModels/incompressible',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/phaseIncompressible',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/twoPhaseInter/incompressibleInterPhaseTransportModel',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_incompressibleInterPhaseTransportModels = library(
    'incompressibleInterPhaseTransportModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_incompressibleInterPhaseTransportModels)
