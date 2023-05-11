srcfiles = [
    lnInclude_hack,
    'decompositionInformation.C',
    'decompositionModel.C',
    'dimFieldDecomposer.C',
    'fvFieldDecomposer.C',
    'fvFieldDecomposerCache.C',
    'pointFieldDecomposer.C',
    'pointFieldDecomposerCache.C',
    'lagrangianFieldDecomposer.C',
    'lagrangianFieldDecomposerCache.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_decompositionMethods,
    lib_lagrangian,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/parallel/decompose/decompose',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/lagrangian/basic',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompose',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_decompose = library(
    'decompose',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_decompose)
