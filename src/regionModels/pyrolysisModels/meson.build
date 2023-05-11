srcfiles = [
    lnInclude_hack,
    'pyrolysisModel/pyrolysisModel.C',
    'pyrolysisModel/pyrolysisModelNew.C',
    'reactingOneDim/reactingOneDim.C',
    'noPyrolysis/noPyrolysis.C',
    'thermo/thermo.C',
    'pyrolysisModel/pyrolysisModelCollection.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_chemistryModel,
    lib_specie,
    lib_fluidThermophysicalModels,
    lib_solidChemistryModel,
    lib_solidThermo,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_regionModels,
    lib_radiationModels,
    lib_reactionThermophysicalModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/regionModels/pyrolysisModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/chemistryModel',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidThermo',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidChemistryModel',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/reactionThermo',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/radiation',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/regionModels/regionModel',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/regionModels/pyrolysisModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_pyrolysisModels = library(
    'pyrolysisModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_pyrolysisModels)
