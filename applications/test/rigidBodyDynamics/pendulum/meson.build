srcfiles = [
    lnInclude_hack,
    'pendulum.C',
]
link_with = [
    lib_OpenFOAM,
    lib_rigidBodyDynamics,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/rigidBodyDynamics/pendulum',
    '-I' + recursive_include_dirs / 'src/rigidBodyDynamics',
    '-I' + recursive_include_dirs / 'applications/test/rigidBodyDynamics/pendulum',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_pendulum = executable(
    'Test_pendulum',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
