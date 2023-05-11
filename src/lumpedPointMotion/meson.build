srcfiles = [
    lnInclude_hack,
    'controller/lumpedPointController.C',
    'state/lumpedPointState.C',
    'state/lumpedPointStateWriter.C',
    'movement/lumpedPointMovement.C',
    'movement/lumpedPointMovementWriter.C',
    'movement/lumpedPointIOMovement.C',
    'tools/lumpedPointTools.C',
    'pointPatchFields/lumpedPointDisplacementPointPatchVectorField.C',
]
link_with = [
    lib_finiteVolume,
    lib_fileFormats,
    lib_meshTools,
    lib_dynamicMesh,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/lumpedPointMotion',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/lumpedPointMotion',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_lumpedPointMotion = library(
    'lumpedPointMotion',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_lumpedPointMotion)
