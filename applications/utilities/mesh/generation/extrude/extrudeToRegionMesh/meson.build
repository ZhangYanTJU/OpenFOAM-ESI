srcfiles = [
    lnInclude_hack,
    'extrudeToRegionMesh.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_dynamicMesh,
    lib_reconstruct,
    lib_extrudeModel,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/mesh/generation/extrude/extrudeToRegionMesh',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/parallel/reconstruct/reconstruct',
    '-I' + recursive_include_dirs / 'src/mesh/extrudeModel',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/generation/extrude/extrudeToRegionMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_extrudeToRegionMesh = executable(
    'extrudeToRegionMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
