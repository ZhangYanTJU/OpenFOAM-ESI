srcfiles = [
    lnInclude_hack,
    flexgen.process('fluent3DMeshToFoam.L'),
]
link_with = [
    lib_OpenFOAM,
    lib_meshTools,
    lib_dynamicMesh,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/mesh/conversion/fluent3DMeshToFoam',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/conversion/fluent3DMeshToFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_fluent3DMeshToFoam = executable(
    'fluent3DMeshToFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
