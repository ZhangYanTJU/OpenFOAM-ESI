srcfiles = [
    lnInclude_hack,
    'twoPhaseEulerFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_sampling,
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_specie,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
    lib_incompressibleTransportModels,
    lib_compressibleTwoPhaseSystem,
    lib_phaseCompressibleTurbulenceModels,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/twoPhaseEulerFoam',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/twoPhaseEuler/phaseCompressibleTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/twoPhaseEuler/twoPhaseSystem',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/transportModels/incompressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/phaseCompressible',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/twoPhaseEulerFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_twoPhaseEulerFoam = executable(
    'twoPhaseEulerFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
