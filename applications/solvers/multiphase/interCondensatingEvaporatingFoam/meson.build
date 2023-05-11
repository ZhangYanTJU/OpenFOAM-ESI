subdir('temperaturePhaseChangeTwoPhaseMixtures')

srcfiles = [
    lnInclude_hack,
    'interCondensatingEvaporatingFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_sampling,
    lib_dynamicFvMesh,
    lib_phaseTemperatureChangeTwoPhaseMixtures,
    lib_twoPhaseMixture,
    lib_twoPhaseProperties,
    lib_interfaceProperties,
    lib_incompressibleTransportModels,
    lib_turbulenceModels,
    lib_incompressibleTurbulenceModels,
    lib_fluidThermophysicalModels,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/interCondensatingEvaporatingFoam',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/interPhaseChangeFoam',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/interFoam',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/VoF',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/interCondensatingEvaporatingFoam/temperaturePhaseChangeTwoPhaseMixtures',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fvOptions',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/transportModels/incompressible',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/interCondensatingEvaporatingFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_interCondensatingEvaporatingFoam = executable(
    'interCondensatingEvaporatingFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
