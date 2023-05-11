srcfiles = [
    lnInclude_hack,
    'helpType/helpType.C',
    'helpType/helpTypeNew.C',
    'helpBoundary/helpBoundary.C',
    'helpFunctionObject/helpFunctionObject.C',
    'helpSolver/helpSolver.C',
    'doxygenXmlParser/doxygenXmlParser.C',
]
link_with = [
    lib_finiteVolume,
    lib_turbulenceModels,
    lib_incompressibleTurbulenceModels,
    lib_compressibleTurbulenceModels,
    lib_thermoTools,
    lib_radiationModels,
    lib_fluidThermophysicalModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/miscellaneous/foamHelp/helpTypes',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'applications/utilities/miscellaneous/foamHelp/helpTypes',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_helpTypes = library(
    'helpTypes',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_helpTypes)
