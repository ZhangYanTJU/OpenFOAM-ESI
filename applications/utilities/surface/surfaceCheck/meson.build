srcfiles = [
    lnInclude_hack,
    'surfaceCheck.C',
]
link_with = [
    lib_OpenFOAM,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/surface/surfaceCheck',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/utilities/surface/surfaceCheck',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_surfaceCheck = executable(
    'surfaceCheck',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
