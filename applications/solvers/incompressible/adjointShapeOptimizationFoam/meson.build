srcfiles = [
    lnInclude_hack,
    'adjointOutletPressure/adjointOutletPressureFvPatchScalarField.C',
    'adjointOutletVelocity/adjointOutletVelocityFvPatchVectorField.C',
    'adjointShapeOptimizationFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_turbulenceModels,
    lib_incompressibleTurbulenceModels,
    lib_incompressibleTransportModels,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/incompressible/adjointShapeOptimizationFoam',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/incompressible',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + meson.source_root() / 'src/transportModels/incompressible/singlePhaseTransportModel',
    '-I' + recursive_include_dirs / 'applications/solvers/incompressible/adjointShapeOptimizationFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_adjointShapeOptimizationFoam = executable(
    'adjointShapeOptimizationFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
