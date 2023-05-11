srcfiles = [
    lnInclude_hack,
    'regionModel/regionModel.C',
    'singleLayerRegion/singleLayerRegion.C',
    'regionModel1D/regionModel1D.C',
    'derivedFvPatches/mappedVariableThicknessWall/mappedVariableThicknessWallFvPatch.C',
    'regionModelFunctionObject/regionModelFunctionObject/regionModelFunctionObject.C',
    'regionModelFunctionObject/regionModelFunctionObject/regionModelFunctionObjectNew.C',
    'regionModelFunctionObject/regionModelFunctionObject/regionModelFunctionObjectList.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/regionModels/regionModel',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/regionModels/regionModel',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_regionModels = library(
    'regionModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_regionModels)
