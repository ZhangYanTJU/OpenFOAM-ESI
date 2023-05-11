subdir('alphaFieldFunctions')

srcfiles = [
    lnInclude_hack,
    'setAlphaField.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_sampling,
    lib_geometricVoF,
    lib_twoPhaseMixture,
    lib_twoPhaseProperties,
    lib_alphaFieldFunctions,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/preProcessing/setAlphaField',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/setAlphaField/alphaFieldFunctions',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/sampling',
    '-I' + recursive_include_dirs / 'src/transportModels/geometricVoF',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/setAlphaField',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_setAlphaField = executable(
    'setAlphaField',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
