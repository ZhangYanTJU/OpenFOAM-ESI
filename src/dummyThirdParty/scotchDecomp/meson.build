srcfiles = [
    lnInclude_hack,
    'dummyScotchDecomp.C',
]
link_with = [
    lib_decompositionMethods,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/dummyThirdParty/scotchDecomp',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/decompositionMethods',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/scotchDecomp',
    '-I' + recursive_include_dirs / 'src/dummyThirdParty/scotchDecomp',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_dummy_slash_scotchDecomp = library(
    'dummy_slash_scotchDecomp',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_dummy_slash_scotchDecomp)
