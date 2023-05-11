srcfiles = [
    lnInclude_hack,
    'DEShybrid/DEShybrid.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_turbulenceModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/TurbulenceModels/schemes',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/schemes',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_turbulenceModelSchemes = library(
    'turbulenceModelSchemes',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_turbulenceModelSchemes)
