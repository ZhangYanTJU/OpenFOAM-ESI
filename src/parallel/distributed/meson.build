srcfiles = [
    lnInclude_hack,
    'distributedTriSurfaceMesh/distributedTriSurfaceMesh.C',
    'patchDistMethods/exact/exactPatchDistMethod.C',
]
link_with = [
    lib_finiteVolume,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_decompositionMethods,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/parallel/distributed',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/parallel/distributed',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_distributed = library(
    'distributed',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_distributed)
