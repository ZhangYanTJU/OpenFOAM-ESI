srcfiles = [
    lnInclude_hack,
    'interfaceTrackingFvMesh.C',
    'freeSurfacePointDisplacement.C',
    'fvPatchFields/freeSurfacePressure/freeSurfacePressureFvPatchScalarField.C',
    'fvPatchFields/freeSurfaceVelocity/freeSurfaceVelocityFvPatchVectorField.C',
    'functionObjects/pointHistory/pointHistory.C',
    'functionObjects/writeFreeSurface/writeFreeSurface.C',
]
link_with = [
    lib_turbulenceModels,
    lib_incompressibleTurbulenceModels,
    lib_incompressibleTransportModels,
    lib_finiteVolume,
    lib_finiteArea,
    lib_fileFormats,
    lib_meshTools,
    lib_dynamicFvMesh,
    lib_fvMotionSolvers,
    lib_dynamicMesh,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/dynamicFaMesh/interfaceTrackingFvMesh',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + meson.source_root() / 'src/transportModels/incompressible/singlePhaseTransportModel',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/finiteArea',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicFvMesh',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/fvMotionSolver',
    '-I' + recursive_include_dirs / 'src/dynamicFaMesh/interfaceTrackingFvMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_interfaceTrackingFvMesh = library(
    'interfaceTrackingFvMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_interfaceTrackingFvMesh)
