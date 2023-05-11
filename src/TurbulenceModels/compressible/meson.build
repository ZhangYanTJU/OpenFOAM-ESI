srcfiles = [
    lnInclude_hack,
    'compressibleTurbulenceModel.C',
    'turbulentFluidThermoModels/turbulentFluidThermoModels.C',
]
link_with = [
    lib_compressibleTransportModels,
    lib_finiteVolume,
    lib_meshTools,
    lib_radiationModels,
    lib_fluidThermophysicalModels,
    lib_solidThermo,
    lib_solidSpecie,
    lib_turbulenceModels,
    lib_specie,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/radiation',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidThermo',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidSpecie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/thermophysicalProperties',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_compressibleTurbulenceModels = library(
    'compressibleTurbulenceModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_compressibleTurbulenceModels)
