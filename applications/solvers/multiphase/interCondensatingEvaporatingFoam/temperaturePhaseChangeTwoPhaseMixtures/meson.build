srcfiles = [
    lnInclude_hack,
    'temperaturePhaseChangeTwoPhaseMixtures/temperaturePhaseChangeTwoPhaseMixture.C',
    'temperaturePhaseChangeTwoPhaseMixtures/temperaturePhaseChangeTwoPhaseMixtureNew.C',
    'thermoIncompressibleTwoPhaseMixture/thermoIncompressibleTwoPhaseMixture.C',
    'twoPhaseMixtureEThermo/twoPhaseMixtureEThermo.C',
    'constant/constant.C',
    'interfaceHeatResistance/interfaceHeatResistance.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_geometricVoF,
    lib_twoPhaseMixture,
    lib_twoPhaseProperties,
    lib_interfaceProperties,
    lib_incompressibleTransportModels,
    lib_fluidThermophysicalModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/interCondensatingEvaporatingFoam/temperaturePhaseChangeTwoPhaseMixtures',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/transportModels/geometricVoF',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/transportModels/incompressible',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/interCondensatingEvaporatingFoam/temperaturePhaseChangeTwoPhaseMixtures',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_phaseTemperatureChangeTwoPhaseMixtures = library(
    'phaseTemperatureChangeTwoPhaseMixtures',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_phaseTemperatureChangeTwoPhaseMixtures)
