srcfiles = [
    lnInclude_hack,
    'topoChangerFvMesh/topoChangerFvMesh.C',
    'rawTopoChangerFvMesh/rawTopoChangerFvMesh.C',
    'dynamicMotionSolverTopoFvMesh/dynamicMotionSolverTopoFvMesh.C',
    'movingConeTopoFvMesh/movingConeTopoFvMesh.C',
]
link_with = [
    lib_finiteVolume,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_dynamicMesh,
    lib_dynamicFvMesh,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/topoChangerFvMesh',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/topoChangerFvMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_topoChangerFvMesh = library(
    'topoChangerFvMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_topoChangerFvMesh)
