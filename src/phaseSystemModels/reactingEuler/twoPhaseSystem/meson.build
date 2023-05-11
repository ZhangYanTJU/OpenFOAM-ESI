srcfiles = [
    lnInclude_hack,
    'twoPhaseSystem.C',
    'twoPhaseSystemNew.C',
    'twoPhaseSystems.C',
    'diameterModels/IATE/IATE.C',
    'diameterModels/IATE/IATEsources/IATEsource/IATEsource.C',
    'diameterModels/IATE/IATEsources/wakeEntrainmentCoalescence/wakeEntrainmentCoalescence.C',
    'diameterModels/IATE/IATEsources/turbulentBreakUp/turbulentBreakUp.C',
    'diameterModels/IATE/IATEsources/randomCoalescence/randomCoalescence.C',
    'diameterModels/IATE/IATEsources/phaseChange/phaseChange.C',
    'diameterModels/IATE/IATEsources/wallBoiling/wallBoiling.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_combustionModels,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_reactingMultiphaseSystem,
    lib_saturationModel,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/phaseSystemModels/reactingEuler/twoPhaseSystem',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/reactingEuler/multiphaseSystem',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/reactingEuler/saturationModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/reactionThermo',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/phaseCompressible',
    '-I' + recursive_include_dirs / 'src/combustionModels',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/reactingEuler/twoPhaseSystem',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_reactingTwoPhaseSystem = library(
    'reactingTwoPhaseSystem',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_reactingTwoPhaseSystem)
