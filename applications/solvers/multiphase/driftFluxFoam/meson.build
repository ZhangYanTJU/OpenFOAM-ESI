subdir('mixtureViscosityModels')
subdir('relativeVelocityModels')

srcfiles = [
    lnInclude_hack,
    'incompressibleTwoPhaseInteractingMixture/incompressibleTwoPhaseInteractingMixture.C',
    'compressibleTurbulenceModels.C',
    'driftFluxFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_sampling,
    lib_compressibleTransportModels,
    lib_incompressibleTransportModels,
    lib_twoPhaseMixture,
    lib_driftFluxTransportModels,
    lib_driftFluxRelativeVelocityModels,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/driftFluxFoam',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/driftFluxFoam/incompressibleTwoPhaseInteractingMixture',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/driftFluxFoam/mixtureViscosityModels',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/driftFluxFoam/relativeVelocityModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/driftFluxFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_driftFluxFoam = executable(
    'driftFluxFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
