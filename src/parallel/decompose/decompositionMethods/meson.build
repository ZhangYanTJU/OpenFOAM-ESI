srcfiles = [
    lnInclude_hack,
    'decompositionMethod/decompositionMethod.C',
    'geomDecomp/geomDecomp.C',
    'simpleGeomDecomp/simpleGeomDecomp.C',
    'hierarchGeomDecomp/hierarchGeomDecomp.C',
    'manualDecomp/manualDecomp.C',
    'multiLevelDecomp/multiLevelDecomp.C',
    'metisLikeDecomp/metisLikeDecomp.C',
    'structuredDecomp/structuredDecomp.C',
    'randomDecomp/randomDecomp.C',
    'noDecomp/noDecomp.C',
    'decompositionConstraints/decompositionConstraint/decompositionConstraint.C',
    'decompositionConstraints/preserveBaffles/preserveBafflesConstraint.C',
    'decompositionConstraints/preserveFaceZones/preserveFaceZonesConstraint.C',
    'decompositionConstraints/preservePatches/preservePatchesConstraint.C',
    'decompositionConstraints/geometric/geometricConstraint.C',
    'decompositionConstraints/singleProcessorFaceSets/singleProcessorFaceSetsConstraint.C',
    'decompositionConstraints/refinementHistory/refinementHistoryConstraint.C',
]
link_with = [
    lib_finiteVolume,
    lib_fileFormats,
    lib_meshTools,
    lib_dynamicMesh,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_decompositionMethods = library(
    'decompositionMethods',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_decompositionMethods)
