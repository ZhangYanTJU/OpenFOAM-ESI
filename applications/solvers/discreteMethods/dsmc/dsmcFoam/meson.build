srcfiles = [
    lnInclude_hack,
    'dsmcFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_meshTools,
    lib_finiteVolume,
    lib_lagrangian,
    lib_DSMC,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/discreteMethods/dsmc/dsmcFoam',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/lagrangian/DSMC',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/solvers/discreteMethods/dsmc/dsmcFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_dsmcFoam = executable(
    'dsmcFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
