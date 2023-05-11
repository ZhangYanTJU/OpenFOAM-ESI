srcfiles = [
    lnInclude_hack,
    'decomposePar.C',
    'domainDecomposition.C',
    'domainDecompositionMesh.C',
    'domainDecompositionDistribute.C',
    'domainDecompositionWrite.C',
    'domainDecompositionDryRun.C',
    'domainDecompositionDryRunWrite.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteArea,
    lib_finiteVolume,
    lib_fileFormats,
    lib_meshTools,
    lib_lagrangian,
    lib_dynamicMesh,
    lib_genericPatchFields,
    lib_decompositionMethods,
    lib_decompose,
    lib_faDecompose,
    lib_kahipDecomp,
    lib_metisDecomp,
    lib_scotchDecomp,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/parallelProcessing/decomposePar',
    '-I' + recursive_include_dirs / 'src/finiteArea',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompose',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/faDecompose',
    '-I' + recursive_include_dirs / 'applications/utilities/parallelProcessing/decomposePar',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_decomposePar = executable(
    'decomposePar',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
