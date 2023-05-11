srcfiles = [
    lnInclude_hack,
    'foamToTetDualMesh.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/postProcessing/dataConversion/foamToTetDualMesh',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'applications/utilities/postProcessing/dataConversion/foamToTetDualMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_foamToTetDualMesh = executable(
    'foamToTetDualMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
