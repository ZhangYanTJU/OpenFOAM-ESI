srcfiles = [
    lnInclude_hack,
    'phaseModel/phaseModel.C',
    'alphaContactAngle/alphaContactAngleFvPatchScalarField.C',
    'multiphaseMixtureThermo.C',
]
link_with = [
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_specie,
    lib_finiteVolume,
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/compressibleMultiphaseInterFoam/multiphaseMixtureThermo',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/compressibleMultiphaseInterFoam/multiphaseMixtureThermo',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_multiphaseMixtureThermo = library(
    'multiphaseMixtureThermo',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_multiphaseMixtureThermo)
