srcfiles = [
    lnInclude_hack,
    'energySpectrum/energySpectrum.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_randomProcesses,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/functionObjects/randomProcesses',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/randomProcesses',
    '-I' + recursive_include_dirs / 'src/functionObjects/randomProcesses',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_randomProcessesFunctionObjects = library(
    'randomProcessesFunctionObjects',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_randomProcessesFunctionObjects)
