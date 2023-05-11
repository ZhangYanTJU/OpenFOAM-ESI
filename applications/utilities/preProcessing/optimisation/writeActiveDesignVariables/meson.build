srcfiles = [
    lnInclude_hack,
    'writeActiveDesignVariables.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_meshTools,
    lib_fvOptions,
    lib_dynamicMesh,
    lib_fvMotionSolvers,
    lib_adjointOptimisation,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/preProcessing/optimisation/writeActiveDesignVariables',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + meson.source_root() / 'src/finiteVolume/cfdTools',
    '-I' + recursive_include_dirs / 'src/dynamicMesh',
    '-I' + recursive_include_dirs / 'src/optimisation/adjointOptimisation/adjoint',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/optimisation/writeActiveDesignVariables',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_writeActiveDesignVariables = executable(
    'writeActiveDesignVariables',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
