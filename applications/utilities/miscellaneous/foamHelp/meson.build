subdir('helpTypes')

srcfiles = [
    lnInclude_hack,
    'foamHelp.C',
]
link_with = [
    lib_OpenFOAM,
    lib_helpTypes,
    lib_finiteVolume,
    lib_meshTools,
    lib_turbulenceModels,
    lib_incompressibleTurbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
    lib_radiationModels,
    lib_fluidThermophysicalModels,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/miscellaneous/foamHelp',
    '-I' + recursive_include_dirs / 'applications/utilities/miscellaneous/foamHelp/helpTypes',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/utilities/miscellaneous/foamHelp',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_foamHelp = executable(
    'foamHelp',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
