srcfiles = [
    lnInclude_hack,
    'makeFaMesh.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteArea,
    lib_fileFormats,
    lib_surfMesh,
    lib_meshTools,
    lib_faDecompose,
    lib_faReconstruct,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/finiteArea/makeFaMesh',
    '-I' + recursive_include_dirs / 'src/finiteArea',
    '-I' + recursive_include_dirs / 'src/fileFormats',
    '-I' + recursive_include_dirs / 'src/surfMesh',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/parallel/decompose/faDecompose',
    '-I' + recursive_include_dirs / 'src/parallel/reconstruct/faReconstruct',
    '-I' + recursive_include_dirs / 'applications/utilities/finiteArea/makeFaMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_makeFaMesh = executable(
    'makeFaMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
