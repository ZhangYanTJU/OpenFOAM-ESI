srcfiles = [
    lnInclude_hack,
    'sizeDistribution/sizeDistribution.C',
    'phaseForces/phaseForces.C',
]
link_with = [
    lib_finiteVolume,
    lib_fieldFunctionObjects,
    lib_reactingMultiphaseSystem,
    lib_dynamicMesh,
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/functionObjects/phaseSystems',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/reactingEuler/multiphaseSystem',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/functionObjects/field',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/functionObjects/phaseSystems',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_phaseFunctionObjects = library(
    'phaseFunctionObjects',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_phaseFunctionObjects)
