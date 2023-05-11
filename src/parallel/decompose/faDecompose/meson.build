srcfiles = [
    lnInclude_hack,
    'faFieldDecomposer.C',
    'faFieldDecomposerCache.C',
    'faMeshDecomposition.C',
]
link_with = [
    lib_finiteArea,
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/parallel/decompose/faDecompose',
    '-I' + recursive_include_dirs / 'src/finiteArea',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/faDecompose',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_faDecompose = library(
    'faDecompose',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_faDecompose)
