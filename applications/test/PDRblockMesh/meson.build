srcfiles = [
    lnInclude_hack,
    'Test-PDRblockMesh.C',
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
    '-DFULLDEBUG',
    '-I' + meson.source_root() / 'applications/test/PDRblockMesh',
    '-I' + recursive_include_dirs / 'src/mesh/blockMesh',
    '-I' + recursive_include_dirs / 'applications/test/PDRblockMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_PDRblockMesh = executable(
    'Test_PDRblockMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
