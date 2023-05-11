srcfiles = [
    lnInclude_hack,
    'mirrorFvMesh.C',
    'mirrorMesh.C',
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
    '-I' + meson.source_root() / 'applications/utilities/mesh/manipulation/mirrorMesh',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/manipulation/mirrorMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_mirrorMesh = executable(
    'mirrorMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
