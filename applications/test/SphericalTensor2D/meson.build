srcfiles = [
    lnInclude_hack,
    'Test-SphericalTensor2D.C',
]
link_with = [
    lib_OpenFOAM,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/SphericalTensor2D',
    '-I' + recursive_include_dirs / 'applications/test/SphericalTensor2D',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_SphericalTensor2D = executable(
    'Test_SphericalTensor2D',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
