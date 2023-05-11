srcfiles = [
    lnInclude_hack,
    'common/ccmBase.C',
    'common/ccmInternal.C',
    'reader/ccmReader.C',
    'reader/ccmReaderAux.C',
    'reader/ccmReaderMesh.C',
    'reader/ccmReaderSolution.C',
    'reader/ccmReaderOptions.C',
    'writer/ccmWriter.C',
    'writer/ccmWriterMesh.C',
    'writer/ccmWriterSolution.C',
]
link_with = [
    lib_finiteVolume,
    lib_fileFormats,
    lib_conversion,
]
dependencies = [
    ccmio_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'src/conversion/ccm',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/conversion',
    '-I' + meson.source_root() / 'src/conversion/ccm',
    '-I' + recursive_include_dirs / 'src/conversion/ccm',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_ccm = library(
    'ccm',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_ccm)
