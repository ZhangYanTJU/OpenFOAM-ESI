srcfiles = [
    lnInclude_hack,
    'common/reader/meshReader.C',
    'common/reader/meshReaderAux.C',
    'common/reader/calcPointCells.C',
    'common/reader/createPolyCells.C',
    'common/reader/createPolyBoundary.C',
    'common/writer/meshWriter.C',
    'common/tables/boundaryRegion.C',
    'common/tables/cellTable.C',
    'fire/FIREMeshReader.C',
    'fire/FIREMeshWriter.C',
    'fire/checkFireEdges.C',
    'starcd/STARCDMeshReader.C',
    'starcd/STARCDMeshWriter.C',
    'polyDualMesh/polyDualMesh.C',
    'vtk/output/foamVtkSurfaceFieldWriter.C',
]
link_with = [
    lib_finiteVolume,
    lib_finiteArea,
    lib_fileFormats,
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/conversion',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/finiteArea',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/conversion',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_conversion = library(
    'conversion',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_conversion)

subdir('ccm')
