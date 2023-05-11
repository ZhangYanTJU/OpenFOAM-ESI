srcfiles = [
    lnInclude_hack,
    'Test-spline.C',
]
link_with = [
    lib_OpenFOAM,
    lib_blockMesh,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/spline',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/mesh/blockMesh',
    '-I' + recursive_include_dirs / 'applications/test/spline',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_spline = executable(
    'Test_spline',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
