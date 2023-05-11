srcfiles = [
    lnInclude_hack,
    'laminarFlameSpeed/laminarFlameSpeed.C',
    'laminarFlameSpeed/laminarFlameSpeedNew.C',
    'constant/constant.C',
    'Gulders/Gulders.C',
    'GuldersEGR/GuldersEGR.C',
    'RaviPetersen/RaviPetersen.C',
]
link_with = [
    lib_finiteVolume,
    lib_compressibleTransportModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/thermophysicalModels/laminarFlameSpeed',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/transportModels/compressible',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/basic',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/reactionThermo',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/laminarFlameSpeed',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_laminarFlameSpeedModels = library(
    'laminarFlameSpeedModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_laminarFlameSpeedModels)
