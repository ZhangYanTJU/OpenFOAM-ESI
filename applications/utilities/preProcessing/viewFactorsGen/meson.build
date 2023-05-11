srcfiles = [
    lnInclude_hack,
    'viewFactorsGen.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
    lib_surfMesh,
    lib_meshTools,
    lib_distributed,
    lib_radiationModels,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-Wno-old-style-cast',
    '-DCGAL_HEADER_ONLY',
    '-I' + meson.source_root() / 'applications/utilities/preProcessing/viewFactorsGen',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/parallel/distributed',
    '-I' + recursive_include_dirs / 'applications/utilities/preProcessing/viewFactorsGen',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
if cgal_dep.found()
    dependencies += cgal_dep
else
    cpp_args += '-DNO_CGAL'
endif
exe_viewFactorsGen = executable(
    'viewFactorsGen',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
