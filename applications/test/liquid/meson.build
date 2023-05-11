srcfiles = [
    lnInclude_hack,
    'Test-liquid.C',
]
link_with = [
    lib_OpenFOAM,
    lib_thermophysicalProperties,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/liquid',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/thermophysicalProperties',
    '-I' + recursive_include_dirs / 'applications/test/liquid',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_liquid = executable(
    'Test_liquid',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
