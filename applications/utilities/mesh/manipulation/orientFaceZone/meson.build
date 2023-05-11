srcfiles = [
    lnInclude_hack,
    'orientFaceZone.C',
]
link_with = [
    lib_OpenFOAM,
    lib_meshTools,
    lib_snappyHexMesh,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/mesh/manipulation/orientFaceZone',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/mesh/snappyHexMesh',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/manipulation/orientFaceZone',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_orientFaceZone = executable(
    'orientFaceZone',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
