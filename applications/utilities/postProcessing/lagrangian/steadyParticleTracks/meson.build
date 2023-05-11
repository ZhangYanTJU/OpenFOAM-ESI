srcfiles = [
    lnInclude_hack,
    'steadyParticleTracks.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fileFormats,
    lib_meshTools,
    lib_lagrangian,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/postProcessing/lagrangian/steadyParticleTracks',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'applications/utilities/postProcessing/lagrangian/steadyParticleTracks',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_steadyParticleTracks = executable(
    'steadyParticleTracks',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
