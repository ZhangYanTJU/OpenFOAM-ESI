srcfiles = [
    lnInclude_hack,
    'reconstructPar.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteArea,
    lib_finiteVolume,
    lib_meshTools,
    lib_lagrangian,
    lib_genericPatchFields,
    lib_dynamicMesh,
    lib_reconstruct,
    lib_faReconstruct,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/parallelProcessing/reconstructPar',
    '-I' + recursive_include_dirs / 'src/finiteArea',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/parallel/reconstruct/reconstruct',
    '-I' + recursive_include_dirs / 'src/parallel/reconstruct/faReconstruct',
    '-I' + recursive_include_dirs / 'applications/utilities/parallelProcessing/reconstructPar',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_reconstructPar = executable(
    'reconstructPar',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
