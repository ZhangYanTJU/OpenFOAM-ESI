srcfiles = [
    lnInclude_hack,
    'reactingTwoPhaseEulerFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_sampling,
    lib_reactingMultiphaseSystem,
    lib_reactingTwoPhaseSystem,
    lib_twoPhaseReactingTurbulenceModels,
    lib_thermoTools,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/reactingTwoPhaseEulerFoam',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/reactingEuler/twoPhaseSystem',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/reactingEuler/twoPhaseCompressibleTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/reactingEuler/multiphaseSystem',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/phaseCompressible',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/reactingTwoPhaseEulerFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_reactingTwoPhaseEulerFoam = executable(
    'reactingTwoPhaseEulerFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
