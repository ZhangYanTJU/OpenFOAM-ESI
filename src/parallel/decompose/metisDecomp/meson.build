srcfiles = [
    lnInclude_hack,
    'metisDecomp.C',
]
link_with = [
    lib_decompositionMethods,
]
dependencies = [
    metis_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'src/parallel/decompose/metisDecomp',
    '-I' + meson.source_root() / 'src/parallel/decompose/metisDecomp',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/metisDecomp',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_metisDecomp = library(
    'metisDecomp',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_metisDecomp)
