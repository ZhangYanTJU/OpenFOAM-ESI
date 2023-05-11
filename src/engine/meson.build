srcfiles = [
    lnInclude_hack,
    'engineTime/engineTime/engineTime.C',
    'engineTime/engineTime/engineTimeNew.C',
    'engineTime/freePiston/freePiston.C',
    'engineTime/crankConRod/crankConRod.C',
    'ignition/ignition.C',
    'ignition/ignitionIO.C',
    'ignition/ignitionSite.C',
    'ignition/ignitionSiteIO.C',
    'engineValve/engineValve.C',
    'enginePiston/enginePiston.C',
    'engineMesh/engineMesh/engineMesh.C',
    'engineMesh/engineMesh/engineMeshNew.C',
    'engineMesh/staticEngineMesh/staticEngineMesh.C',
    'engineMesh/layeredEngineMesh/layeredEngineMesh.C',
    'engineMesh/fvMotionSolverEngineMesh/fvMotionSolverEngineMesh.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_dynamicMesh,
    lib_fvMotionSolvers,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/engine',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/fvMotionSolver',
    '-I' + recursive_include_dirs / 'src/engine',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_engine = library(
    'engine',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_engine)
