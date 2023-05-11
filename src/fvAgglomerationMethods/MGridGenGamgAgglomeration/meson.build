srcfiles = [
    lnInclude_hack,
    'MGridGenGAMGAgglomeration.C',
    'MGridGenGAMGAgglomerate.C',
]
link_with = []
dependencies = [
    mgrid_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'src/fvAgglomerationMethods/MGridGenGamgAgglomeration',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + meson.source_root() / 'src/fvAgglomerationMethods/MGridGenGamgAgglomeration',
    '-I' + recursive_include_dirs / 'src/fvAgglomerationMethods/MGridGenGamgAgglomeration',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_MGridGenGAMGAgglomeration = library(
    'MGridGenGAMGAgglomeration',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_MGridGenGAMGAgglomeration)
