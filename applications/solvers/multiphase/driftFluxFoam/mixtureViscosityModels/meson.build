srcfiles = [
    lnInclude_hack,
    'mixtureViscosityModel/mixtureViscosityModel.C',
    'mixtureViscosityModel/mixtureViscosityModelNew.C',
    'plastic/plastic.C',
    'BinghamPlastic/BinghamPlastic.C',
    'slurry/slurry.C',
]
link_with = [
    lib_twoPhaseMixture,
    lib_incompressibleTransportModels,
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/driftFluxFoam/mixtureViscosityModels',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/driftFluxFoam/incompressibleTwoPhaseInteractingMixture',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/transportModels/incompressible',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/driftFluxFoam/mixtureViscosityModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_driftFluxTransportModels = library(
    'driftFluxTransportModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_driftFluxTransportModels)
