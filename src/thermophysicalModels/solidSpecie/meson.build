srcfiles = [
    lnInclude_hack,
    'reaction/reactions/makeSolidReactions.C',
]
link_with = [
    lib_specie,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/thermophysicalModels/solidSpecie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/thermophysicalProperties',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidSpecie',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_solidSpecie = library(
    'solidSpecie',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_solidSpecie)
