srcfiles = [
    lnInclude_hack,
    'foamCalc.C',
]
link_with = [
    lib_OpenFOAM,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/tools/foamCalc',
    '-I' + recursive_include_dirs / 'applications/tools/foamCalc',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_foamCalc = executable(
    'foamCalc',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
