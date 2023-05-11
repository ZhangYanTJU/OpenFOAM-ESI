srcfiles = [
    lnInclude_hack,
    'chemistryModel/basicChemistryModel/basicChemistryModel.C',
    'chemistryModel/BasicChemistryModel/BasicChemistryModels.C',
    'chemistryModel/TDACChemistryModel/reduction/makeChemistryReductionMethods.C',
    'chemistryModel/TDACChemistryModel/tabulation/makeChemistryTabulationMethods.C',
    'chemistrySolver/chemistrySolver/makeChemistrySolvers.C',
    'functionObjects/specieReactionRates/specieReactionRates.C',
    'functionObjects/BilgerMixtureFraction/BilgerMixtureFraction.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_ODE,
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_reactionThermophysicalModels,
    lib_specie,
    lib_thermophysicalProperties,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/thermophysicalModels/chemistryModel',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/ODE',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/reactionThermo',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/thermophysicalProperties',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/chemistryModel',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_chemistryModel = library(
    'chemistryModel',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_chemistryModel)
