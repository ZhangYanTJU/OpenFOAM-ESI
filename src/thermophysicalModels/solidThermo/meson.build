srcfiles = [
    lnInclude_hack,
    'solidThermo/solidThermo.C',
    'solidThermo/solidThermos.C',
    'solidReactionThermo/solidReactionThermo.C',
    'solidReactionThermo/solidReactionThermos.C',
    'solidThermo/heZoneSolidThermos.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_specie,
    lib_solidSpecie,
    lib_fluidThermophysicalModels,
    lib_reactionThermophysicalModels,
    lib_thermophysicalProperties,
    lib_compressibleTransportModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/thermophysicalModels/solidThermo',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidSpecie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/reactionThermo',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/thermophysicalProperties',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidThermo',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_solidThermo = library(
    'solidThermo',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_solidThermo)
