srcfiles = [
    lnInclude_hack,
    'forces/forces.C',
    'forceCoeffs/forceCoeffs.C',
    'propellerInfo/propellerInfo.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_surfMesh,
    lib_fluidThermophysicalModels,
    lib_incompressibleTransportModels,
    lib_compressibleTransportModels,
    lib_turbulenceModels,
    lib_incompressibleTurbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_specie,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/functionObjects/forces',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'src/functionObjects/forces',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_forces = library(
    'forces',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_forces)
