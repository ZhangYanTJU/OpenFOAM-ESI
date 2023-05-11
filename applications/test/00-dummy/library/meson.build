subdir('dummy')

srcfiles = [
    lnInclude_hack,
    'dummyLib.cpp',
]
link_with = []
dependencies = []
cpp_args = [
    '-DWM_ARCH=""',
    '-DWM_COMPILER=""',
    '-DWM_COMPILE_OPTION=""',
    '-DWM_OPTIONS=""',
    '-I' + meson.source_root() / 'applications/test/00-dummy/library',
    '-I' + recursive_include_dirs / 'applications/test/00-dummy/library',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_TestDummy_OpenFOAM = library(
    'TestDummy_OpenFOAM',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
pkg.generate(lib_TestDummy_OpenFOAM)

subdir('mpi')
