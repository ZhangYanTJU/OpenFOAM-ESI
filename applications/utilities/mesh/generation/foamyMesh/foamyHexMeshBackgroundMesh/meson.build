srcfiles = [
    lnInclude_hack,
    'foamyHexMeshBackgroundMesh.C',
]
link_with = [
    lib_OpenFOAM,
    lib_conformalVoronoiMesh,
    lib_decompositionMethods,
    lib_decompose,
    lib_meshTools,
    lib_dynamicMesh,
    lib_sampling,
    lib_finiteVolume,
]
dependencies = [
    m_dep,
    dl_dep,
    cgal_dep,
    mpfr_dep,
    gmp_dep,
]
cpp_args = [
    '-DNDEBUG',
    '-I' + meson.source_root() / 'applications/utilities/mesh/generation/foamyMesh/foamyHexMeshBackgroundMesh',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/generation/foamyMesh/conformalVoronoiMesh',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompose',
    '-I' + recursive_include_dirs / 'src/mesh/snappyHexMesh',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/generation/foamyMesh/foamyHexMeshBackgroundMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_foamyHexMeshBackgroundMesh = executable(
    'foamyHexMeshBackgroundMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
