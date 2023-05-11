srcfiles = [
    lnInclude_hack,
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
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/Pstream/dummy',
    '-I' + recursive_include_dirs / 'src/Pstream/dummy',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_dummy_slash_Pstream = library(
    'dummy_slash_Pstream',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_dummy_slash_Pstream)
