srcfiles = [
    lnInclude_hack,
    'Test-fieldMapping.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_dynamicMesh,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-DFULLDEBUG',
    '-g',
    '-O0',
    '-I' + meson.source_root() / 'applications/test/fieldMapping',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/test/fieldMapping',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_fieldMapping = executable(
    'Test_fieldMapping',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
