srcfiles = [
    lnInclude_hack,
    'setSet.C',
]
link_with = [
    lib_OpenFOAM,
    lib_meshTools,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/mesh/manipulation/setSet',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/manipulation/setSet',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
if readline_dep.found()
    cpp_args += '-DHAVE_LIBREADLINE'
    dependencies += readline_dep
endif
exe_setSet = executable(
    'setSet',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
