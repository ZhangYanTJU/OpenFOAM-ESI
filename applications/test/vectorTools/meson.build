srcfiles = [
    lnInclude_hack,
    'Test-vectorTools.C',
]
link_with = [
    lib_OpenFOAM,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/vectorTools',
    '-I' + meson.source_root() / 'applications/utilities/mesh/generation/foamyMesh/conformalVoronoiMesh/vectorTools',
    '-I' + recursive_include_dirs / 'applications/test/vectorTools',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_vectorTools = executable(
    'Test_vectorTools',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
