srcfiles = [
    lnInclude_hack,
    'overBuoyantPimpleDyMFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_sampling,
    lib_meshTools,
    lib_fvOptions,
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_radiationModels,
    lib_specie,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
    lib_atmosphericModels,
    lib_overset,
    lib_dynamicFvMesh,
    lib_topoChangerFvMesh,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/heatTransfer/buoyantPimpleFoam/overBuoyantPimpleDyMFoam',
    '-I' + meson.source_root() / 'applications/solvers/heatTransfer/buoyantPimpleFoam',
    '-I' + meson.source_root() / 'applications/solvers/compressible/rhoPimpleFoam/overRhoPimpleDyMFoam',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/overset',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/radiation',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'applications/solvers/heatTransfer/buoyantPimpleFoam/overBuoyantPimpleDyMFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_overBuoyantPimpleDyMFoam = executable(
    'overBuoyantPimpleDyMFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
