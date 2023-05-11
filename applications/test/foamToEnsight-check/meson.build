srcfiles = [
    lnInclude_hack,
    'foamToEnsight-check.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_finiteArea,
    lib_fileFormats,
    lib_meshTools,
    lib_conversion,
    lib_lagrangianIntermediate,
    lib_genericPatchFields,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/foamToEnsight-check',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/finiteArea',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/conversion',
    '-I' + recursive_include_dirs / 'src/lagrangian/intermediate',
    '-I' + recursive_include_dirs / 'applications/test/foamToEnsight-check',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_foamToEnsight_check = executable(
    'foamToEnsight_check',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
