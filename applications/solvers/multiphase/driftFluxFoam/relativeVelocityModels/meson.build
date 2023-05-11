srcfiles = [
    lnInclude_hack,
    'relativeVelocityModel/relativeVelocityModel.C',
    'simple/simple.C',
    'general/general.C',
]
link_with = [
    lib_incompressibleTransportModels,
    lib_finiteVolume,
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/driftFluxFoam/relativeVelocityModels',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/driftFluxFoam/incompressibleTwoPhaseInteractingMixture',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/driftFluxFoam/mixtureViscosityModels',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/driftFluxFoam/relativeVelocityModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_driftFluxRelativeVelocityModels = library(
    'driftFluxRelativeVelocityModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_driftFluxRelativeVelocityModels)
