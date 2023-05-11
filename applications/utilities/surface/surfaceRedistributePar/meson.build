srcfiles = [
    lnInclude_hack,
    'surfaceRedistributePar.C',
]
link_with = [
    lib_OpenFOAM,
    lib_surfMesh,
    lib_meshTools,
    lib_distributed,
    lib_decompose,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/surface/surfaceRedistributePar',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompose',
    '-I' + recursive_include_dirs / 'src/parallel/distributed',
    '-I' + recursive_include_dirs / 'applications/utilities/surface/surfaceRedistributePar',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_surfaceRedistributePar = executable(
    'surfaceRedistributePar',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
