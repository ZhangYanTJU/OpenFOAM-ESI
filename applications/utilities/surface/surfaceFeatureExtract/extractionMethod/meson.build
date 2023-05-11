srcfiles = [
    lnInclude_hack,
    'surfaceFeaturesExtraction.C',
    'extractFromFile.C',
    'extractFromNone.C',
    'extractFromSurface.C',
]
link_with = [
    lib_surfMesh,
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/surface/surfaceFeatureExtract/extractionMethod',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/utilities/surface/surfaceFeatureExtract/extractionMethod',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_surfaceFeatureExtract = library(
    'surfaceFeatureExtract',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_surfaceFeatureExtract)
