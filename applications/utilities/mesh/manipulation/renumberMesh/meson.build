srcfiles = [
    lnInclude_hack,
    'renumberMesh.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_meshTools,
    lib_dynamicMesh,
    lib_genericPatchFields,
    lib_renumberMethods,
    lib_reconstruct,
    lib_decompositionMethods,
    lib_kahipDecomp,
    lib_metisDecomp,
    lib_scotchDecomp,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/mesh/manipulation/renumberMesh',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/renumber/renumberMethods',
    '-I' + recursive_include_dirs / 'src/renumber/zoltanRenumber',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/parallel/reconstruct/reconstruct',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/manipulation/renumberMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
if zoltan_dep.found()
    cpp_args += '-DHAVE_ZOLTAN'
    dependencies += zoltan_dep
endif
exe_renumberMesh = executable(
    'renumberMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
