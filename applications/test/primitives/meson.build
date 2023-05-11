srcfiles = [
    lnInclude_hack,
    'Test-primitives.C',
]
link_with = [
    lib_OpenFOAM,
    lib_fileFormats,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-DCOMPAT_OPENFOAM_ORG',
    '-I' + meson.source_root() / 'applications/test/primitives',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'applications/test/primitives',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_primitives = executable(
    'Test_primitives',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
