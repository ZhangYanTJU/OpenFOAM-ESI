srcfiles = [
    lnInclude_hack,
    'Test-fft.C',
]
link_with = [
    lib_OpenFOAM,
    lib_randomProcesses,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/fft',
    '-I' + recursive_include_dirs / 'src/randomProcesses',
    '-I' + recursive_include_dirs / 'applications/test/fft',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_fft = executable(
    'Test_fft',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
