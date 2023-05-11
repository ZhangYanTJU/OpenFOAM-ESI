srcfiles = [
    lnInclude_hack,
    'foamUpgradeCyclics.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_genericPatchFields,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/preProcessing/foamUpgradeCyclics',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/foamUpgradeCyclics',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_foamUpgradeCyclics = executable(
    'foamUpgradeCyclics',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
