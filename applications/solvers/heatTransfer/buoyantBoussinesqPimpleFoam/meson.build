srcfiles = [
    lnInclude_hack,
    'buoyantBoussinesqPimpleFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_dynamicMesh,
    lib_dynamicFvMesh,
    lib_sampling,
    lib_turbulenceModels,
    lib_incompressibleTurbulenceModels,
    lib_incompressibleTransportModels,
    lib_radiationModels,
    lib_atmosphericModels,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/heatTransfer/buoyantBoussinesqPimpleFoam',
    '-I' + meson.source_root() / 'applications/solvers/heatTransfer/buoyantBoussinesqSimpleFoam',
    '-I' + meson.source_root() / 'applications/solvers/incompressible/pimpleFoam',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + meson.source_root() / 'src/transportModels/incompressible/singlePhaseTransportModel',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/radiation',
    '-I' + recursive_include_dirs / 'applications/solvers/heatTransfer/buoyantBoussinesqPimpleFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_buoyantBoussinesqPimpleFoam = executable(
    'buoyantBoussinesqPimpleFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
