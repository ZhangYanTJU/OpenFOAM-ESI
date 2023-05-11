srcfiles = [
    lnInclude_hack,
    'zoltanRenumber.C',
]
link_with = [
    lib_meshTools,
]
dependencies = [
    zoltan_dep,
]
cpp_args = [
    '-DMPICH_SKIP_MPICXX',
    '-DOMPI_SKIP_MPICXX',
    '-I' + meson.source_root() / 'src/renumber/zoltanRenumber',
    '-I' + meson.source_root() / 'src/renumber/zoltanRenumber',
    '-I' + recursive_include_dirs / 'src/renumber/renumberMethods',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/renumber/zoltanRenumber',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_zoltanRenumber = library(
    'zoltanRenumber',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_zoltanRenumber)
