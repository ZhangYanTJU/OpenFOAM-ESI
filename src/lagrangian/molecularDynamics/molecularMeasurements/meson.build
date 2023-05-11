srcfiles = [
    lnInclude_hack,
    'distribution/distribution.C',
]
link_with = []
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/lagrangian/molecularDynamics/molecularMeasurements',
    '-I' + recursive_include_dirs / 'src/lagrangian/molecularDynamics/molecularMeasurements',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_molecularMeasurements = library(
    'molecularMeasurements',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_molecularMeasurements)
