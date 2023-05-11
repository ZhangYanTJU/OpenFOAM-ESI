srcfiles = [
    lnInclude_hack,
    'immiscibleIncompressibleTwoPhaseMixture.C',
]
link_with = [
    lib_twoPhaseMixture,
    lib_twoPhaseProperties,
    lib_incompressibleTransportModels,
    lib_interfaceProperties,
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/transportModels/immiscibleIncompressibleTwoPhaseMixture',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/incompressible',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/transportModels/immiscibleIncompressibleTwoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_immiscibleIncompressibleTwoPhaseMixture = library(
    'immiscibleIncompressibleTwoPhaseMixture',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_immiscibleIncompressibleTwoPhaseMixture)
