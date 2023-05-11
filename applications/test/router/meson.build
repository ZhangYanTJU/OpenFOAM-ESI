srcfiles = [
    lnInclude_hack,
    'router.C',
    'Test-processorRouter.C',
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
    '-I' + meson.source_root() / 'applications/test/router',
    '-I' + meson.source_root() / 'applications/test/router/Gather',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'applications/test/router',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_processorRouter = executable(
    'Test_processorRouter',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
