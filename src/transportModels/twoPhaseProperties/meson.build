srcfiles = [
    lnInclude_hack,
    'alphaContactAngle/constantAlphaContactAngle/constantAlphaContactAngleFvPatchScalarField.C',
    'alphaContactAngle/temperatureDependentAlphaContactAngle/temperatureDependentAlphaContactAngleFvPatchScalarField.C',
    'alphaContactAngle/dynamicAlphaContactAngle/dynamicAlphaContactAngleFvPatchScalarField.C',
    'alphaContactAngle/timeVaryingAlphaContactAngle/timeVaryingAlphaContactAngleFvPatchScalarField.C',
    'alphaFixedPressure/alphaFixedPressureFvPatchScalarField.C',
]
link_with = [
    lib_twoPhaseMixture,
    lib_interfaceProperties,
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/transportModels/twoPhaseProperties',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseProperties',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_twoPhaseProperties = library(
    'twoPhaseProperties',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_twoPhaseProperties)
