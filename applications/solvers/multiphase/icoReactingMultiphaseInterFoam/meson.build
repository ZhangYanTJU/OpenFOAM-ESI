subdir('laserDTRM')

srcfiles = [
    lnInclude_hack,
    'icoReactingMultiphaseInterFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_sampling,
    lib_dynamicFvMesh,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
    lib_compressibleTransportModels,
    lib_radiationModels,
    lib_fluidThermophysicalModels,
    lib_solidThermo,
    lib_solidSpecie,
    lib_twoPhaseProperties,
    lib_laserDTRM,
    lib_incompressibleMultiphaseSystems,
    lib_compressibleMultiPhaseTurbulenceModels,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/icoReactingMultiphaseInterFoam',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/multiphaseInter/phasesSystem',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/multiphaseInter/compressibleMultiPhaseTurbulenceModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fvOptions',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/radiation',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/reactionThermo',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/interFoam',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/icoReactingMultiphaseInterFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_icoReactingMultiphaseInterFoam = executable(
    'icoReactingMultiphaseInterFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
