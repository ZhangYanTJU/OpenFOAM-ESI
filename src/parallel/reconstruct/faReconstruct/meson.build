srcfiles = [
    lnInclude_hack,
    'processorFaMeshes.C',
    'faFieldReconstructor.C',
    'faMeshReconstructor.C',
]
link_with = [
    lib_finiteArea,
    lib_finiteVolume,
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/parallel/reconstruct/faReconstruct',
    '-I' + recursive_include_dirs / 'src/finiteArea',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/parallel/reconstruct/faReconstruct',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_faReconstruct = library(
    'faReconstruct',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_faReconstruct)
