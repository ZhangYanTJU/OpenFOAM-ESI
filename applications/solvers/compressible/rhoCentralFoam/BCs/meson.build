srcfiles = [
    lnInclude_hack,
    'U/maxwellSlipUFvPatchVectorField.C',
    'T/smoluchowskiJumpTFvPatchScalarField.C',
    'rho/fixedRhoFvPatchScalarField.C',
]
link_with = [
    lib_finiteVolume,
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_specie,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/compressible/rhoCentralFoam/BCs',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'applications/solvers/compressible/rhoCentralFoam/BCs',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_rhoCentralFoam = library(
    'rhoCentralFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_rhoCentralFoam)
