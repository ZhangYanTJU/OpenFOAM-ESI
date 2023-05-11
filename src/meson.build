subdir('ODE')
subdir('OSspecific')
subdir('Pstream')

srcfiles = [
    lnInclude_hack,
    'lagrangian/distributionModels/distributionModel/distributionModel.C',
    'lagrangian/distributionModels/distributionModel/distributionModelNew.C',
    'lagrangian/distributionModels/binned/binned.C',
    'lagrangian/distributionModels/exponential/exponential.C',
    'lagrangian/distributionModels/fixedValue/fixedValue.C',
    'lagrangian/distributionModels/general/general.C',
    'lagrangian/distributionModels/multiNormal/multiNormal.C',
    'lagrangian/distributionModels/normal/normal.C',
    'lagrangian/distributionModels/RosinRammler/RosinRammler.C',
    'lagrangian/distributionModels/massRosinRammler/massRosinRammler.C',
    'lagrangian/distributionModels/uniform/uniform.C',
]
link_with = []
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/lagrangian/distributionModels',
    '-I' + recursive_include_dirs / 'src/lagrangian/distributionModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_distributionModels = library(
    'distributionModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_distributionModels)

subdir('rigidBodyDynamics')
subdir('OpenFOAM')
subdir('fileFormats')
subdir('surfMesh')
subdir('meshTools')
subdir('finiteArea')
subdir('finiteVolume')

srcfiles = [
    lnInclude_hack,
    'lagrangian/basic/particle/particle.C',
    'lagrangian/basic/particle/particleIO.C',
    'lagrangian/basic/passiveParticle/passiveParticleCloud.C',
    'lagrangian/basic/indexedParticle/indexedParticleCloud.C',
    'lagrangian/basic/injectedParticle/injectedParticle.C',
    'lagrangian/basic/injectedParticle/injectedParticleIO.C',
    'lagrangian/basic/injectedParticle/injectedParticleCloud.C',
    'lagrangian/basic/InteractionLists/referredWallFace/referredWallFace.C',
]
link_with = [
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_lagrangian = library(
    'lagrangian',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_lagrangian)

subdir('mesh')
subdir('conversion')
subdir('fvAgglomerationMethods')
subdir('genericPatchFields')
subdir('sampling')
subdir('dynamicMesh')
subdir('randomProcesses')
subdir('dynamicFvMesh')
subdir('lumpedPointMotion')
subdir('parallel')
subdir('transportModels')
subdir('topoChangerFvMesh')
subdir('dummyThirdParty')
subdir('overset')
subdir('renumber')
subdir('thermophysicalModels')
subdir('TurbulenceModels')
subdir('combustionModels')
subdir('fvOptions')
subdir('phaseSystemModels')
subdir('atmosphericModels')
subdir('fvMotionSolver')
subdir('waveModels')
subdir('thermoTools')
subdir('dynamicFaMesh')
subdir('engine')

srcfiles = [
    lnInclude_hack,
    'mesh/snappyHexMesh/snappyHexMeshDriver/snappyLayerDriver.C',
    'mesh/snappyHexMesh/snappyHexMeshDriver/snappyLayerDriverSinglePass.C',
    'mesh/snappyHexMesh/snappyHexMeshDriver/snappySnapDriver.C',
    'mesh/snappyHexMesh/snappyHexMeshDriver/snappySnapDriverFeature.C',
    'mesh/snappyHexMesh/snappyHexMeshDriver/snappyRefineDriver.C',
    'mesh/snappyHexMesh/snappyHexMeshDriver/snappyVoxelMeshDriver.C',
    'mesh/snappyHexMesh/snappyHexMeshDriver/layerParameters/layerParameters.C',
    'mesh/snappyHexMesh/snappyHexMeshDriver/refinementParameters/refinementParameters.C',
    'mesh/snappyHexMesh/snappyHexMeshDriver/snapParameters/snapParameters.C',
    'mesh/snappyHexMesh/meshRefinement/meshRefinementBaffles.C',
    'mesh/snappyHexMesh/meshRefinement/meshRefinement.C',
    'mesh/snappyHexMesh/meshRefinement/meshRefinementMerge.C',
    'mesh/snappyHexMesh/meshRefinement/meshRefinementProblemCells.C',
    'mesh/snappyHexMesh/meshRefinement/meshRefinementRefine.C',
    'mesh/snappyHexMesh/meshRefinement/meshRefinementGapRefine.C',
    'mesh/snappyHexMesh/meshRefinement/meshRefinementBlock.C',
    'mesh/snappyHexMesh/meshRefinement/wallPoints.C',
    'mesh/snappyHexMesh/meshRefinement/patchFaceOrientation.C',
    'mesh/snappyHexMesh/meshRefinement/weightedPosition.C',
    'mesh/snappyHexMesh/refinementFeatures/refinementFeatures.C',
    'mesh/snappyHexMesh/refinementSurfaces/surfaceZonesInfo.C',
    'mesh/snappyHexMesh/refinementSurfaces/refinementSurfaces.C',
    'mesh/snappyHexMesh/shellSurfaces/shellSurfaces.C',
    'mesh/snappyHexMesh/trackedParticle/trackedParticle.C',
    'mesh/snappyHexMesh/trackedParticle/trackedParticleCloud.C',
    'mesh/snappyHexMesh/externalDisplacementMeshMover/displacementMeshMoverMotionSolver.C',
    'mesh/snappyHexMesh/externalDisplacementMeshMover/externalDisplacementMeshMover.C',
    'mesh/snappyHexMesh/externalDisplacementMeshMover/medialAxisMeshMover.C',
    'mesh/snappyHexMesh/externalDisplacementMeshMover/displacementMotionSolverMeshMover.C',
    'mesh/snappyHexMesh/externalDisplacementMeshMover/zeroFixedValue/zeroFixedValuePointPatchFields.C',
    'mesh/snappyHexMesh/externalDisplacementMeshMover/fieldSmoother/fieldSmoother.C',
]
link_with = [
    lib_finiteVolume,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_dynamicMesh,
    lib_sampling,
    lib_lagrangian,
    lib_fvMotionSolvers,
    lib_overset,
    lib_distributed,
    lib_decompositionMethods,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/mesh/snappyHexMesh',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/overset',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/parallel/distributed',
    '-I' + recursive_include_dirs / 'src/mesh/snappyHexMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_snappyHexMesh = library(
    'snappyHexMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_snappyHexMesh)

subdir('optimisation')
subdir('faOptions')
subdir('regionModels')
subdir('regionFaModels')
subdir('lagrangian')
subdir('functionObjects')
subdir('rigidBodyMeshMotion')
subdir('sixDoFRigidBodyMotion')
subdir('sixDoFRigidBodyState')
