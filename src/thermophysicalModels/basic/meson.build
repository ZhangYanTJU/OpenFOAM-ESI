srcfiles = [
    lnInclude_hack,
    'basicThermo/basicThermo.C',
    'fluidThermo/fluidThermo.C',
    'psiThermo/psiThermo.C',
    'psiThermo/psiThermos.C',
    'psiThermo/psiZoneThermos.C',
    'rhoThermo/rhoThermo.C',
    'rhoThermo/rhoThermos.C',
    'rhoThermo/liquidThermo.C',
    'rhoThermo/rhoZoneThermos.C',
    'derivedFvPatchFields/fixedEnergy/fixedEnergyFvPatchScalarField.C',
    'derivedFvPatchFields/gradientEnergy/gradientEnergyFvPatchScalarField.C',
    'derivedFvPatchFields/mixedEnergy/mixedEnergyFvPatchScalarField.C',
    'derivedFvPatchFields/energyJump/energyJump/energyJumpFvPatchScalarField.C',
    'derivedFvPatchFields/energyJump/energyJumpAMI/energyJumpAMIFvPatchScalarField.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_compressibleTransportModels,
    lib_specie,
    lib_thermophysicalProperties,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/thermophysicalProperties',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_fluidThermophysicalModels = library(
    'fluidThermophysicalModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_fluidThermophysicalModels)
