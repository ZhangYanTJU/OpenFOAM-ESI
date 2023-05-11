srcfiles = [
    lnInclude_hack,
    'makeWave.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fileFormats,
    lib_meshTools,
    lib_dynamicMesh,
    lib_dynamicFvMesh,
    lib_fvMotionSolvers,
    lib_interfaceTrackingFvMesh,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'tutorials/incompressible/pimpleFoam/laminar/sloshing2D/makeWave',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/fvMotionSolver',
    '-I' + recursive_include_dirs / 'src/dynamicFaMesh/interfaceTrackingFvMesh',
    '-I' + recursive_include_dirs / 'tutorials/incompressible/pimpleFoam/laminar/sloshing2D/makeWave',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_makeWave = executable(
    'makeWave',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
