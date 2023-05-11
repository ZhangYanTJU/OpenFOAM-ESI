srcfiles = [
    lnInclude_hack,
    'computeSensitivities.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_sampling,
    lib_turbulenceModels,
    lib_incompressibleTransportModels,
    lib_incompressibleTurbulenceModels,
    lib_fileFormats,
    lib_surfMesh,
    lib_adjointOptimisation,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/postProcessing/optimisation/computeSensitivities',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fvOptions',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + meson.source_root() / 'src/transportModels/incompressible/singlePhaseTransportModel',
    '-I' + recursive_include_dirs / 'src/optimisation/adjointOptimisation/adjoint',
    '-I' + recursive_include_dirs / 'applications/utilities/postProcessing/optimisation/computeSensitivities',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_computeSensitivities = executable(
    'computeSensitivities',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
