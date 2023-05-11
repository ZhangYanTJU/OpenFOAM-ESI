subdir('extractionMethod')

srcfiles = [
    lnInclude_hack,
    'surfaceFeatureExtract.C',
]
link_with = [
    lib_OpenFOAM,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_surfaceFeatureExtract,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/surface/surfaceFeatureExtract',
    '-I' + recursive_include_dirs / 'applications/utilities/surface/surfaceFeatureExtract/extractionMethod',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/utilities/surface/surfaceFeatureExtract',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_surfaceFeatureExtract = executable(
    'surfaceFeatureExtract',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
