srcfiles = [
    lnInclude_hack,
    'barotropicCompressibilityModel/barotropicCompressibilityModel.C',
    'barotropicCompressibilityModel/barotropicCompressibilityModelNew.C',
    'linear/linear.C',
    'Wallis/Wallis.C',
    'Chung/Chung.C',
]
link_with = [
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/thermophysicalModels/barotropicCompressibilityModel',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/barotropicCompressibilityModel',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_barotropicCompressibilityModel = library(
    'barotropicCompressibilityModel',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_barotropicCompressibilityModel)
