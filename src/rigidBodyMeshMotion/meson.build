srcfiles = [
    lnInclude_hack,
    'rigidBodyMeshMotion/rigidBodyMeshMotion.C',
    'rigidBodyMeshMotionSolver/rigidBodyMeshMotionSolver.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_rigidBodyDynamics,
    lib_forces,
    lib_dynamicMesh,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/rigidBodyMeshMotion',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/rigidBodyDynamics',
    '-I' + recursive_include_dirs / 'src/functionObjects/forces',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/rigidBodyMeshMotion',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_rigidBodyMeshMotion = library(
    'rigidBodyMeshMotion',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_rigidBodyMeshMotion)
