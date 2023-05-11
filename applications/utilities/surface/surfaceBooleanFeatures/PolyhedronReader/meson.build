srcfiles = [
    lnInclude_hack,
    'PolyhedronReader.C',
]
link_with = [
    lib_surfMesh,
    lib_meshTools,
]
dependencies = [
    cgal_dep,
    mpfr_dep,
    gmp_dep,
]
cpp_args = [
    '-DNDEBUG',
    '-I' + meson.source_root() / 'applications/utilities/surface/surfaceBooleanFeatures/PolyhedronReader',
    '-I' + meson.source_root() / 'applications/utilities/surface/surfaceBooleanFeatures',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/utilities/surface/surfaceBooleanFeatures/PolyhedronReader',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_PolyhedronReader = library(
    'PolyhedronReader',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_PolyhedronReader)
