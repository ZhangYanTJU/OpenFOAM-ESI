srcfiles = [
    lnInclude_hack,
    'Test-PatchFunction1.C',
]
link_with = [
    lib_OpenFOAM,
    lib_lagrangianIntermediate,
    lib_radiationModels,
    lib_regionModels,
    lib_finiteVolume,
    lib_meshTools,
    lib_sampling,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-DFULLDEBUG',
    '-g',
    '-O0',
    '-I' + meson.source_root() / 'applications/test/PatchFunction1',
    '-I' + recursive_include_dirs / 'src/lagrangian/intermediate',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/radiation',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/test/PatchFunction1',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_Test_PatchFunction1 = executable(
    'Test_PatchFunction1',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
