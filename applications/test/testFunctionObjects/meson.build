srcfiles = [
    lnInclude_hack,
    'fakeError/fakeErrorFunctionObject.C',
]
link_with = []
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/test/testFunctionObjects',
    '-I' + recursive_include_dirs / 'applications/test/testFunctionObjects',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_testFunctionObjects = library(
    'testFunctionObjects',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
pkg.generate(lib_testFunctionObjects)
