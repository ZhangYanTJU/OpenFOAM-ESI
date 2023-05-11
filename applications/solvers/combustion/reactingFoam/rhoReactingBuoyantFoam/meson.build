srcfiles = [
    lnInclude_hack,
    'rhoReactingBuoyantFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_sampling,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
    lib_reactionThermophysicalModels,
    lib_specie,
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_chemistryModel,
    lib_ODE,
    lib_combustionModels,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/combustion/reactingFoam/rhoReactingBuoyantFoam',
    '-I' + meson.source_root() / 'applications/solvers/combustion/reactingFoam',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/reactionThermo',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/chemistryModel',
    '-I' + recursive_include_dirs / 'src/ODE',
    '-I' + recursive_include_dirs / 'src/combustionModels',
    '-I' + recursive_include_dirs / 'applications/solvers/combustion/reactingFoam/rhoReactingBuoyantFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_rhoReactingBuoyantFoam = executable(
    'rhoReactingBuoyantFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
