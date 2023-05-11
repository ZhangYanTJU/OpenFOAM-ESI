srcfiles = [
    lnInclude_hack,
    'foamCellZoneToVTK.C',
]
link_with = [
    lib_OpenFOAM,
    lib_fileFormats,
    lib_meshTools,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/foamCellZoneToVTK',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/test/foamCellZoneToVTK',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_foamCellZoneToVTK = executable(
    'foamCellZoneToVTK',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
