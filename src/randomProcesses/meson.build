srcfiles = [
    lnInclude_hack,
    'Kmesh/Kmesh.C',
    'fft/fft.C',
    'fft/calcEk.C',
    'fft/kShellIntegration.C',
    'processes/UOprocess/UOprocess.C',
    'turbulence/turbGen.C',
    'noise/noiseFFT/noiseFFT.C',
    'noise/noiseModels/noiseModel/noiseModel.C',
    'noise/noiseModels/noiseModel/noiseModelNew.C',
    'noise/noiseModels/pointNoise/pointNoise.C',
    'noise/noiseModels/surfaceNoise/surfaceNoise.C',
    'windowModels/windowModel/windowModel.C',
    'windowModels/windowModel/windowModelNew.C',
    'windowModels/Hanning/Hanning.C',
    'windowModels/uniform/uniform.C',
]
link_with = [
    lib_finiteVolume,
    lib_fileFormats,
    lib_meshTools,
    lib_surfMesh,
    lib_sampling,
]
dependencies = [
    fftw3_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'src/randomProcesses',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + meson.source_root() / 'src/randomProcesses',
    '-I' + recursive_include_dirs / 'src/randomProcesses',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_randomProcesses = library(
    'randomProcesses',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_randomProcesses)
