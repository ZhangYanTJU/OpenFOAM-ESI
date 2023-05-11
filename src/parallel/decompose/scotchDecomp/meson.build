srcfiles = [
    lnInclude_hack,
    'scotchDecomp.C',
]
link_with = [
    lib_decompositionMethods,
]
dependencies = [
    scotch_dep,
    scotcherrexit_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'src/parallel/decompose/scotchDecomp',
    '-I' + meson.source_root() / 'src/parallel/decompose/scotchDecomp',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/scotchDecomp',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_scotchDecomp = library(
    'scotchDecomp',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_scotchDecomp)
