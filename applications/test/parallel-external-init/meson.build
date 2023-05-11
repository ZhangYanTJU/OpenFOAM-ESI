srcfiles = [
    lnInclude_hack,
    'Test-parallel-external-init.C',
]
link_with = [
    lib_OpenFOAM,
]
dependencies = [
    m_dep,
    dl_dep,
    mpi_dep,
]
cpp_args = [
    '-DMPICH_SKIP_MPICXX',
    '-DOMPI_SKIP_MPICXX',
    '-I' + meson.source_root() / 'applications/test/parallel-external-init',
    '-I' + recursive_include_dirs / 'applications/test/parallel-external-init',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_parallel_external_init = executable(
    'Test_parallel_external_init',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
