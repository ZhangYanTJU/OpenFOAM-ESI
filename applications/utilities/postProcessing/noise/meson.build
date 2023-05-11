srcfiles = [
    lnInclude_hack,
    'noise.C',
]
link_with = [
    lib_OpenFOAM,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_sampling,
    lib_randomProcesses,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/postProcessing/noise',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/randomProcesses',
    '-I' + meson.source_root() / 'applications/utilities/postProcessing/noise',
    '-I' + recursive_include_dirs / 'applications/utilities/postProcessing/noise',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_noise = executable(
    'noise',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
