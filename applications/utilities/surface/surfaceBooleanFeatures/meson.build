subdir('PolyhedronReader')

srcfiles = [
    lnInclude_hack,
    'surfaceBooleanFeatures.C',
]
link_with = [
    lib_OpenFOAM,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
]
dependencies = [
    m_dep,
    dl_dep,
    cgal_dep,
    mpfr_dep,
    gmp_dep,
]
cpp_args = [
    '-DNDEBUG',
    '-I' + meson.source_root() / 'applications/utilities/surface/surfaceBooleanFeatures',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'applications/utilities/surface/surfaceBooleanFeatures',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
if cgal_dep.found()
    cpp_args += '-I' + meson.source_root() / 'applications/utilities/surface/surfaceBooleanFeatures/PolyhedronReader'
    link_with += lib_PolyhedronReader
    dependencies += cgal_dep
else
    cpp_args += '-DNO_CGAL'
endif
exe_surfaceBooleanFeatures = executable(
    'surfaceBooleanFeatures',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
