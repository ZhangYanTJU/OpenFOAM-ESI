srcfiles = [
    lnInclude_hack,
    'faceSelection/faceSelection.C',
    'faceSelection/faceZoneSelection.C',
    'faceSelection/searchableSurfaceSelection.C',
    'createBaffles.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_meshTools,
    lib_dynamicMesh,
    lib_reconstruct,
    lib_genericPatchFields,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/mesh/manipulation/createBaffles',
    '-I' + meson.source_root() / 'applications/utilities/mesh/manipulation/createBaffles/faceSelection',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/parallel/reconstruct/reconstruct',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/manipulation/createBaffles',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_createBaffles = executable(
    'createBaffles',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
