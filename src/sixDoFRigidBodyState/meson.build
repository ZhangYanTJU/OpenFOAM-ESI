srcfiles = [
    lnInclude_hack,
    'sixDoFRigidBodyState.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_dynamicFvMesh,
    lib_sixDoFRigidBodyMotion,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/sixDoFRigidBodyState',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/sixDoFRigidBodyMotion',
    '-I' + recursive_include_dirs / 'src/sixDoFRigidBodyState',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_sixDoFRigidBodyState = library(
    'sixDoFRigidBodyState',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_sixDoFRigidBodyState)
