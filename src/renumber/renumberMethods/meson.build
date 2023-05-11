srcfiles = [
    lnInclude_hack,
    'renumberMethod/renumberMethod.C',
    'manualRenumber/manualRenumber.C',
    'CuthillMcKeeRenumber/CuthillMcKeeRenumber.C',
    'randomRenumber/randomRenumber.C',
    'springRenumber/springRenumber.C',
    'structuredRenumber/structuredRenumber.C',
    'structuredRenumber/OppositeFaceCellWaveBase.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_decompositionMethods,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/renumber/renumberMethods',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/renumber/renumberMethods',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_renumberMethods = library(
    'renumberMethods',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_renumberMethods)
