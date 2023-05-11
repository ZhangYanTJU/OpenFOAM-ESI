subdir('tabulatedWallFunction')

srcfiles = [
    lnInclude_hack,
    'wallFunctionTable.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_meshTools,
    lib_tabulatedWallFunctions,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/preProcessing/wallFunctionTable',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/wallFunctionTable/tabulatedWallFunction',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/wallFunctionTable',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_wallFunctionTable = executable(
    'wallFunctionTable',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
