srcfiles = [
    lnInclude_hack,
    'viscosityModels/viscosityModel/viscosityModel.C',
    'viscosityModels/viscosityModel/viscosityModelNew.C',
    'viscosityModels/Newtonian/Newtonian.C',
    'viscosityModels/powerLaw/powerLaw.C',
    'viscosityModels/CrossPowerLaw/CrossPowerLaw.C',
    'viscosityModels/BirdCarreau/BirdCarreau.C',
    'viscosityModels/HerschelBulkley/HerschelBulkley.C',
    'viscosityModels/Casson/Casson.C',
    'viscosityModels/strainRateFunction/strainRateFunction.C',
    'viscosityModels/Arrhenius/Arrheniuss.C',
    'transportModel/transportModel.C',
    'singlePhaseTransportModel/singlePhaseTransportModel.C',
    'incompressibleTwoPhaseMixture/incompressibleTwoPhaseMixture.C',
]
link_with = [
    lib_twoPhaseMixture,
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/transportModels/incompressible',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/transportModels/incompressible',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_incompressibleTransportModels = library(
    'incompressibleTransportModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_incompressibleTransportModels)
