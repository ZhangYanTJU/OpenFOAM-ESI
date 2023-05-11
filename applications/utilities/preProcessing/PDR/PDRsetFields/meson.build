srcfiles = [
    lnInclude_hack,
    'PDRsetFields.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_blockMesh,
    lib_pdrFields,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/preProcessing/PDR/PDRsetFields',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/PDR/pdrFields',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/mesh/blockMesh',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/PDR/PDRsetFields',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_PDRsetFields = executable(
    'PDRsetFields',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
