srcfiles = [
    lnInclude_hack,
    'processorMeshes.C',
    'fvFieldReconstructor.C',
    'pointFieldReconstructor.C',
    'passivePositionParticleCloud.C',
    'lagrangianReconstructor.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_lagrangian,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/parallel/reconstruct/reconstruct',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/parallel/reconstruct/reconstruct',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_reconstruct = library(
    'reconstruct',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_reconstruct)
