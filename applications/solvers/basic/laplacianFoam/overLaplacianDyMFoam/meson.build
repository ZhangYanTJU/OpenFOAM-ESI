srcfiles = [
    lnInclude_hack,
    'overLaplacianDyMFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_meshTools,
    lib_dynamicFvMesh,
    lib_overset,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/basic/laplacianFoam/overLaplacianDyMFoam',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/overset',
    '-I' + recursive_include_dirs / 'src/overset/include',
    '-I' + recursive_include_dirs / 'applications/solvers/basic/laplacianFoam/overLaplacianDyMFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_overLaplacianDyMFoam = executable(
    'overLaplacianDyMFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
