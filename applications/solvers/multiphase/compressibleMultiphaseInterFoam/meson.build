subdir('multiphaseMixtureThermo')

srcfiles = [
    lnInclude_hack,
    'compressibleMultiphaseInterFoam.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fvOptions,
    lib_meshTools,
    lib_multiphaseMixtureThermo,
    lib_compressibleTransportModels,
    lib_fluidThermophysicalModels,
    lib_specie,
    lib_interfaceProperties,
    lib_turbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/compressibleMultiphaseInterFoam',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/VoF',
    '-I' + meson.source_root() / 'applications/solvers/multiphase/interFoam',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/compressibleMultiphaseInterFoam/multiphaseMixtureThermo',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/transportModels/twoPhaseMixture',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/turbulenceModels',
    '-I' + recursive_include_dirs / 'src/TurbulenceModels/compressible',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/compressibleMultiphaseInterFoam',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_compressibleMultiphaseInterFoam = executable(
    'compressibleMultiphaseInterFoam',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
