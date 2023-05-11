srcfiles = [
    lnInclude_hack,
    'patchSummary.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_meshTools,
    lib_genericPatchFields,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/miscellaneous/patchSummary',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/utilities/miscellaneous/patchSummary',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_patchSummary = executable(
    'patchSummary',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
