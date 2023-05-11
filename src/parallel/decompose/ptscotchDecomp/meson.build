srcfiles = [
    lnInclude_hack,
    'ptscotchDecomp.C',
]
link_with = [
    lib_decompositionMethods,
]
dependencies = [
    scotch_dep,
    ptscotch_dep,
    ptscotcherrexit_dep,
]
cpp_args = [
    '-DMPICH_SKIP_MPICXX',
    '-DOMPI_SKIP_MPICXX',
    '-I' + meson.source_root() / 'src/parallel/decompose/ptscotchDecomp',
    '-I' + meson.source_root() / 'src/parallel/decompose/ptscotchDecomp',
    '-I' + meson.source_root() / 'src/parallel/decompose/ptscotchDecomp',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/ptscotchDecomp',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_ptscotchDecomp = library(
    'ptscotchDecomp',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_ptscotchDecomp)
