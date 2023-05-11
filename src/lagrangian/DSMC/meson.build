srcfiles = [
    lnInclude_hack,
    'parcels/derived/dsmcParcel/dsmcParcel.C',
    'clouds/baseClasses/DSMCBaseCloud/DSMCBaseCloud.C',
    'parcels/derived/dsmcParcel/defineDSMCParcel.C',
    'parcels/derived/dsmcParcel/makeDSMCParcelBinaryCollisionModels.C',
    'parcels/derived/dsmcParcel/makeDSMCParcelWallInteractionModels.C',
    'parcels/derived/dsmcParcel/makeDSMCParcelInflowBoundaryModels.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_lagrangian,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/lagrangian/DSMC',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/lagrangian/DSMC',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_DSMC = library(
    'DSMC',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_DSMC)
