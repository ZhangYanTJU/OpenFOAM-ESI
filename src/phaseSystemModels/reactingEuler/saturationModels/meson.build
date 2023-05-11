srcfiles = [
    lnInclude_hack,
    'saturationModel/saturationModel.C',
    'Antoine/Antoine.C',
    'AntoineExtended/AntoineExtended.C',
    'ArdenBuck/ArdenBuck.C',
    'polynomial/polynomial.C',
    'function1/function1.C',
    'constantSaturationConditions/constantSaturationConditions.C',
]
link_with = [
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/phaseSystemModels/reactingEuler/saturationModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/reactingEuler/saturationModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_saturationModel = library(
    'saturationModel',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_saturationModel)
