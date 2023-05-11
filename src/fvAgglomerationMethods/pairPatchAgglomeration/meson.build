srcfiles = [
    lnInclude_hack,
    'pairPatchAgglomeration.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_meshTools,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/fvAgglomerationMethods/pairPatchAgglomeration',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/fvAgglomerationMethods/pairPatchAgglomeration',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_pairPatchAgglomeration = library(
    'pairPatchAgglomeration',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_pairPatchAgglomeration)
