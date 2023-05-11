srcfiles = [
    lnInclude_hack,
    'SloanRenumber.C',
]
link_with = [
    lib_meshTools,
    lib_decompositionMethods,
    lib_renumberMethods,
]
dependencies = [
    boost_system_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'src/renumber/SloanRenumber',
    '-I' + meson.source_root() / 'src/renumber/SloanRenumber',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/renumber/renumberMethods',
    '-I' + recursive_include_dirs / 'src/renumber/SloanRenumber',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_SloanRenumber = library(
    'SloanRenumber',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_SloanRenumber)
