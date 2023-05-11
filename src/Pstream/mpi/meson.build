srcfiles = [
    lnInclude_hack,
    'PstreamGlobals.C',
    'UPstream.C',
    'UPstreamAllToAll.C',
    'UPstreamBroadcast.C',
    'UPstreamGatherScatter.C',
    'UPstreamReduce.C',
    'UIPstreamRead.C',
    'UOPstreamWrite.C',
    'UIPBstreamRead.C',
    'UOPBstreamWrite.C',
]
link_with = []
dependencies = [
    mpi_dep,
]
cpp_args = [
    '-DMPICH_SKIP_MPICXX',
    '-DOMPI_SKIP_MPICXX',
    '-I' + meson.source_root() / 'src/Pstream/mpi',
    '-I' + recursive_include_dirs / 'src/Pstream/mpi',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_Pstream = library(
    'Pstream',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_Pstream)
