srcfiles = [
    lnInclude_hack,
    'interfaceProperties.C',
    'interfaceCompression/interfaceCompression.C',
    'alphaContactAngle/alphaContactAngleTwoPhaseFvPatchScalarField.C',
    'surfaceTensionModels/surfaceTensionModel/surfaceTensionModel.C',
    'surfaceTensionModels/surfaceTensionModel/surfaceTensionModelNew.C',
    'surfaceTensionModels/constant/constantSurfaceTension.C',
    'surfaceTensionModels/temperatureDependent/temperatureDependentSurfaceTension.C',
]
link_with = [
    lib_twoPhaseMixture,
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_interfaceProperties = library(
    'interfaceProperties',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_interfaceProperties)
