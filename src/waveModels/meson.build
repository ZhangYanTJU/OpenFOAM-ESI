srcfiles = [
    lnInclude_hack,
    'waveModel/waveModel.C',
    'waveModel/waveModelNew.C',
    'waveGenerationModels/base/waveGenerationModel/waveGenerationModel.C',
    'waveGenerationModels/base/irregularWaveModel/irregularWaveModel.C',
    'waveGenerationModels/base/regularWaveModel/regularWaveModel.C',
    'waveGenerationModels/base/solitaryWaveModel/solitaryWaveModel.C',
    'waveGenerationModels/derived/Boussinesq/BoussinesqWaveModel.C',
    'waveGenerationModels/derived/cnoidal/cnoidalWaveModel.C',
    'waveGenerationModels/derived/Grimshaw/GrimshawWaveModel.C',
    'waveGenerationModels/derived/McCowan/McCowanWaveModel.C',
    'waveGenerationModels/derived/streamFunction/streamFunctionWaveModel.C',
    'waveGenerationModels/derived/StokesII/StokesIIWaveModel.C',
    'waveGenerationModels/derived/StokesI/StokesIWaveModel.C',
    'waveGenerationModels/derived/StokesV/StokesVWaveModel.C',
    'waveGenerationModels/derived/irregularMultiDirectional/irregularMultiDirectionalWaveModel.C',
    'waveAbsorptionModels/base/waveAbsorptionModel/waveAbsorptionModel.C',
    'waveAbsorptionModels/derived/shallowWaterAbsorption/shallowWaterAbsorption.C',
    'derivedFvPatchFields/waveVelocity/waveVelocityFvPatchVectorField.C',
    'derivedFvPatchFields/waveAlpha/waveAlphaFvPatchScalarField.C',
    'derivedPointPatchFields/waveMaker/waveMakerPointPatchVectorField.C',
    'fvOptions/multiphaseMangrovesSource/multiphaseMangrovesSource.C',
    'fvOptions/multiphaseMangrovesTurbulenceModel/multiphaseMangrovesTurbulenceModel.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_fvOptions,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/waveModels',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/fvOptions',
    '-I' + recursive_include_dirs / 'src/waveModels',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_waveModels = library(
    'waveModels',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_waveModels)
