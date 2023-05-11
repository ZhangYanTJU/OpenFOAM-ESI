srcfiles = [
    lnInclude_hack,
    'PDRarrays.C',
    'PDRarraysAnalyse.C',
    'PDRarraysCalc.C',
    'PDRmeshArrays.C',
    'PDRparams.C',
    'PDRpatchDef.C',
    'PDRlegacyMeshSpec.C',
    'PDRutilsIntersect.C',
    'PDRutilsOverlap.C',
    'obstacles/PDRobstacle.C',
    'obstacles/PDRobstacleIO.C',
    'obstacles/PDRobstacleTypes.C',
    'obstacles/PDRobstacleLegacyIO.C',
    'obstacles/PDRobstacleLegacyRead.C',
]
link_with = [
    lib_finiteVolume,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_blockMesh,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/preProcessing/PDR/pdrFields',
    '-I' + meson.source_root() / 'applications/utilities/preProcessing/PDR/pdrFields/obstacles',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/mesh/blockMesh',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/PDR/pdrFields',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_pdrFields = library(
    'pdrFields',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_pdrFields)
