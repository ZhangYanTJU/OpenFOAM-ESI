srcfiles = [
    lnInclude_hack,
    'phase/phase.C',
    'alphaContactAngle/alphaContactAngleFvPatchScalarField.C',
    'multiphaseMixture.C',
]
link_with = [
    lib_interfaceProperties,
    lib_incompressibleTransportModels,
    lib_finiteVolume,
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/multiphaseInterFoam/multiphaseMixture',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/multiphaseInterFoam/multiphaseMixture/alphaContactAngle',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/incompressible',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/multiphaseInterFoam/multiphaseMixture',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_multiphaseInterFoam = library(
    'multiphaseInterFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_multiphaseInterFoam)
