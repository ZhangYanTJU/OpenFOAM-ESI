srcfiles = [
    lnInclude_hack,
    'surfactantFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteArea,
    lib_finiteVolume,
    lib_meshTools,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/finiteArea/surfactantFoam',
    '-I' + recursive_include_dirs / 'src/finiteArea',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/solvers/finiteArea/surfactantFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_surfactantFoam = executable(
    'surfactantFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
