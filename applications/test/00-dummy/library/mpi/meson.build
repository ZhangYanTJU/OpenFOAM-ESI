srcfiles = [
    lnInclude_hack,
    'dummyMpiLib.cpp',
]
link_with = []
dependencies = [
    mpi_dep,
]
cpp_args = [
    '-DMPICH_SKIP_MPICXX',
    '-DOMPI_SKIP_MPICXX',
    '-DFOAM_MPI=""',
    '-I' + meson.source_root() / 'applications/test/00-dummy/library/mpi',
    '-I' + meson.source_root() / 'applications/test/00-dummy/library',
    '-I' + recursive_include_dirs / 'applications/test/00-dummy/library/mpi',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_TestDummyMpi = library(
    'TestDummyMpi',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
pkg.generate(lib_TestDummyMpi)
