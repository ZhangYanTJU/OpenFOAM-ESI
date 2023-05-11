srcfiles = [
    lnInclude_hack,
    'dynamicFvMesh/dynamicFvMesh.C',
    'dynamicFvMesh/dynamicFvMeshNew.C',
    'staticFvMesh/staticFvMesh.C',
    'dynamicMotionSolverFvMesh/dynamicMotionSolverFvMesh.C',
    'dynamicMultiMotionSolverFvMesh/dynamicMultiMotionSolverFvMesh.C',
    'dynamicInkJetFvMesh/dynamicInkJetFvMesh.C',
    'dynamicRefineFvMesh/dynamicRefineFvMesh.C',
    'dynamicMotionSolverListFvMesh/dynamicMotionSolverListFvMesh.C',
    'simplifiedDynamicFvMesh/simplifiedDynamicFvMeshes.C',
    'simplifiedDynamicFvMesh/simplifiedDynamicFvMesh.C',
    'dynamicMotionSolverFvMeshAMI/dynamicMotionSolverFvMeshAMI.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_dynamicMesh,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_dynamicFvMesh = library(
    'dynamicFvMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_dynamicFvMesh)
