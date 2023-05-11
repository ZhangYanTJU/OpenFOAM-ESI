srcfiles = [
    lnInclude_hack,
    'basicSolidChemistryModel/basicSolidChemistryModel.C',
    'basicSolidChemistryModel/basicSolidChemistryModelNew.C',
    'basicSolidChemistryModel/basicSolidChemistryModels.C',
    'solidChemistrySolver/makeSolidChemistrySolvers.C',
]
link_with = [
    lib_finiteVolume,
    lib_ODE,
    lib_specie,
    lib_chemistryModel,
    lib_fluidThermophysicalModels,
    lib_solidThermo,
    lib_reactionThermophysicalModels,
    lib_compressibleTransportModels,
    lib_thermophysicalProperties,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/thermophysicalModels/solidChemistryModel',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/ODE',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidSpecie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/chemistryModel',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidThermo',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/reactionThermo',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/thermophysicalProperties',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidChemistryModel',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_solidChemistryModel = library(
    'solidChemistryModel',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_solidChemistryModel)
