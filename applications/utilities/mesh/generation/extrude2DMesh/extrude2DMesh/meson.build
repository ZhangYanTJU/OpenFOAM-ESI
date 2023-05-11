srcfiles = [
    lnInclude_hack,
    'extrude2DMesh/extrude2DMesh.C',
    'patchToPoly2DMesh/patchToPoly2DMesh.C',
]
link_with = [
    lib_finiteVolume,
    lib_surfMesh,
    lib_meshTools,
    lib_dynamicMesh,
    lib_extrudeModel,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/mesh/generation/extrude2DMesh/extrude2DMesh',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/mesh/extrudeModel',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/generation/extrude2DMesh/extrude2DMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_extrude2DMesh = library(
    'extrude2DMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_extrude2DMesh)
