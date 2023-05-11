srcfiles = [
    lnInclude_hack,
    'incompressibleTurbulenceModel.C',
    'incompressibleRhoTurbulenceModel.C',
    'turbulentTransportModels/turbulentTransportModels.C',
    'turbulentTransportModels/RAS/qZeta/qZeta.C',
    'turbulentTransportModels/RAS/kkLOmega/kkLOmega.C',
    'turbulentTransportModels/RAS/LamBremhorstKE/LamBremhorstKE.C',
    'turbulentTransportModels/RAS/LienLeschziner/LienLeschziner.C',
    'turbulentTransportModels/RAS/ShihQuadraticKE/ShihQuadraticKE.C',
    'turbulentTransportModels/RAS/LienCubicKE/LienCubicKE.C',
    'turbulentTransportModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatJayatillekeWallFunction/alphatJayatillekeWallFunctionFvPatchScalarField.C',
]
link_with = [
    lib_turbulenceModels,
    lib_finiteVolume,
    lib_meshTools,
    lib_incompressibleTransportModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/TurbulenceModels/incompressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_incompressibleTurbulenceModels = library(
    'incompressibleTurbulenceModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_incompressibleTurbulenceModels)
