srcfiles = [
    lnInclude_hack,
    'Test-fieldDependency.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/fieldDependency',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/test/fieldDependency',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_fieldDependency = executable(
    'Test_fieldDependency',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
