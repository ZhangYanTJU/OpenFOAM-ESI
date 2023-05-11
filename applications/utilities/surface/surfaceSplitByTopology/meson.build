srcfiles = [
    lnInclude_hack,
    'surfaceSplitByTopology.C',
]
link_with = [
    lib_OpenFOAM,
    lib_surfMesh,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/surface/surfaceSplitByTopology',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'applications/utilities/surface/surfaceSplitByTopology',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_surfaceSplitByTopology = executable(
    'surfaceSplitByTopology',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
