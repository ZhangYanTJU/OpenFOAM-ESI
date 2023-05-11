srcfiles = [
    lnInclude_hack,
    'dummyMGridGen.C',
]
link_with = []
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/dummyThirdParty/MGridGen',
    '-I' + recursive_include_dirs / 'src/dummyThirdParty/MGridGen',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_dummy_slash_MGridGen = library(
    'dummy_slash_MGridGen',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_dummy_slash_MGridGen)
