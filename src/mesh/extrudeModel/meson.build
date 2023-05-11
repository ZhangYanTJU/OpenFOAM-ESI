srcfiles = [
    lnInclude_hack,
    'extrudeModel/extrudeModel.C',
    'extrudeModel/extrudeModelNew.C',
    'linearNormal/linearNormal.C',
    'planeExtrusion/planeExtrusion.C',
    'linearDirection/linearDirection.C',
    'linearRadial/linearRadial.C',
    'offsetSurface/offsetSurface.C',
    'radial/radial.C',
    'sigmaRadial/sigmaRadial.C',
    'sector/sector.C',
    'cyclicSector/cyclicSector.C',
    'wedge/wedge.C',
    'polyline/polyline.C',
]
link_with = [
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_blockMesh,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/mesh/extrudeModel',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/mesh/blockMesh',
    '-I' + recursive_include_dirs / 'src/mesh/extrudeModel',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_extrudeModel = library(
    'extrudeModel',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_extrudeModel)
