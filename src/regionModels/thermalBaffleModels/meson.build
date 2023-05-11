srcfiles = [
    lnInclude_hack,
    'thermalBaffleModel/thermalBaffleModel.C',
    'thermalBaffleModel/thermalBaffleModelNew.C',
    'thermalBaffle/thermalBaffle.C',
    'noThermo/noThermo.C',
    'derivedFvPatchFields/thermalBaffle/thermalBaffleFvPatchScalarField.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_compressibleTransportModels,
    lib_regionModels,
    lib_solidThermo,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
    lib_radiationModels,
    lib_dynamicMesh,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/regionModels/thermalBaffleModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/regionModels/regionModel',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/solidThermo',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/radiation',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermoTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/mesh/extrudeModel',
    '-I' + recursive_include_dirs / 'src/regionModels/thermalBaffleModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_thermalBaffleModels = library(
    'thermalBaffleModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_thermalBaffleModels)
