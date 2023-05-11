srcfiles = [
    lnInclude_hack,
    'Test-dictionary2.C',
]
link_with = [
    lib_OpenFOAM,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-DCOMPAT_OPENFOAM_ORG',
    '-I' + meson.source_root() / 'applications/test/dictionary2',
    '-I' + recursive_include_dirs / 'applications/test/dictionary2',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_dictionary2 = executable(
    'Test_dictionary2',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
