srcfiles = [
    lnInclude_hack,
    'laserDTRM.C',
    'DTRMParticle/DTRMParticle.C',
    'DTRMParticle/DTRMParticleIO.C',
    'localDensityAbsorptionEmission/localDensityAbsorptionEmission.C',
    'reflectionModel/reflectionModel/reflectionModel.C',
    'reflectionModel/reflectionModel/reflectionModelNew.C',
    'reflectionModel/noReflection/noReflection.C',
    'reflectionModel/FresnelLaser/FresnelLaser.C',
    'reflectionModel/Fresnel/Fresnel.C',
]
link_with = [
    lib_finiteVolume,
    lib_fileFormats,
    lib_meshTools,
    lib_lagrangian,
    lib_radiationModels,
    lib_incompressibleMultiphaseSystems,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/solvers/multiphase/icoReactingMultiphaseInterFoam/laserDTRM',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/multiphaseInter/phasesSystem',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/radiation',
    '-I' + recursive_include_dirs / 'applications/solvers/multiphase/icoReactingMultiphaseInterFoam/laserDTRM',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_laserDTRM = library(
    'laserDTRM',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_laserDTRM)
