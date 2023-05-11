srcfiles = [
    lnInclude_hack,
    'reducedUnits/reducedUnits.C',
    'reducedUnits/reducedUnitsIO.C',
    'molecule/molecule.C',
    'molecule/moleculeIO.C',
    'moleculeCloud/moleculeCloud.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_lagrangian,
    lib_potential,
    lib_molecularMeasurements,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/lagrangian/molecularDynamics/molecule',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/lagrangian/molecularDynamics/potential',
    '-I' + recursive_include_dirs / 'src/lagrangian/molecularDynamics/molecularMeasurements',
    '-I' + recursive_include_dirs / 'src/lagrangian/molecularDynamics/molecule',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_molecule = library(
    'molecule',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_molecule)
