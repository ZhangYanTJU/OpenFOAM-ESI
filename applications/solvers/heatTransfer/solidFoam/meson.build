srcfiles = [
    lnInclude_hack,
    'solidRegionDiffNo.C',
    'solidFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_sampling,
    lib_dynamicMesh,
    lib_dynamicFvMesh,
    lib_topoChangerFvMesh,
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_solidThermo,
    lib_radiationModels,
    lib_specie,
    lib_thermoTools,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/heatTransfer/solidFoam',
    '-I' + meson.source_root() / 'applications/solvers/heatTransfer/chtMultiRegionFoam/include',
    '-I' + meson.source_root() / 'applications/solvers/heatTransfer/chtMultiRegionFoam/solid',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidThermo',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/radiation',
    '-I' + recursive_include_dirs / 'applications/solvers/heatTransfer/solidFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_solidFoam = executable(
    'solidFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
