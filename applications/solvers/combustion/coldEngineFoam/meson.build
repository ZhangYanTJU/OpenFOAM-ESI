srcfiles = [
    lnInclude_hack,
    'coldEngineFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_engine,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_specie,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/combustion/coldEngineFoam',
    '-I' + meson.source_root() / 'applications/solvers/combustion/XiFoam/XiEngineFoam',
    '-I' + meson.source_root() / 'applications/solvers/combustion/XiFoam',
    '-I' + meson.source_root() / 'applications/solvers/compressible/rhoPimpleFoam',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/engine',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/reactionThermo',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'applications/solvers/combustion/coldEngineFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_coldEngineFoam = executable(
    'coldEngineFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
