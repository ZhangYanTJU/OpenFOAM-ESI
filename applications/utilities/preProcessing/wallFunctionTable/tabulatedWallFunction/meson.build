srcfiles = [
    lnInclude_hack,
    'tabulatedWallFunction/tabulatedWallFunction.C',
    'tabulatedWallFunction/tabulatedWallFunctionNew.C',
    'SpaldingsLaw/SpaldingsLaw.C',
    'general/general.C',
]
link_with = [
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/preProcessing/wallFunctionTable/tabulatedWallFunction',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/wallFunctionTable/tabulatedWallFunction',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_tabulatedWallFunctions = library(
    'tabulatedWallFunctions',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_tabulatedWallFunctions)
