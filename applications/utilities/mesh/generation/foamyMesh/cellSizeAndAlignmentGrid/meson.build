srcfiles = [
    lnInclude_hack,
    'cellSizeAndAlignmentGrid.C',
]
link_with = [
    lib_OpenFOAM,
    lib_conformalVoronoiMesh,
    lib_finiteVolume,
    lib_meshTools,
    lib_decompositionMethods,
    lib_kahipDecomp,
    lib_metisDecomp,
    lib_ptscotchDecomp,
    lib_scotchDecomp,
    lib_dynamicMesh,
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
    '-I' + meson.source_root() / 'applications/utilities/mesh/generation/foamyMesh/cellSizeAndAlignmentGrid',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/mesh/snappyHexMesh',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/generation/foamyMesh/conformalVoronoiMesh',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/generation/foamyMesh/cellSizeAndAlignmentGrid',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_cellSizeAndAlignmentGrid = executable(
    'cellSizeAndAlignmentGrid',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
