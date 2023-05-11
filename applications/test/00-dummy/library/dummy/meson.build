srcfiles = [
    lnInclude_hack,
    'dummyMpiLib.cpp',
]
link_with = []
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/00-dummy/library/dummy',
    '-I' + meson.source_root() / 'applications/test/00-dummy/library',
    '-I' + recursive_include_dirs / 'applications/test/00-dummy/library/dummy',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_dummy_slash_TestDummyMpi = library(
    'dummy_slash_TestDummyMpi',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
pkg.generate(lib_dummy_slash_TestDummyMpi)
