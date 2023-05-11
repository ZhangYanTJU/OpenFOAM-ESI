srcfiles = [
    lnInclude_hack,
    'genericPatchFieldBase/genericPatchFieldBase.C',
    'genericFaPatchField/genericFaPatchFields.C',
    'genericFvPatchField/genericFvPatchFields.C',
    'genericFvsPatchField/genericFvsPatchFields.C',
    'genericPointPatchField/genericPointPatchFields.C',
]
link_with = [
    lib_finiteArea,
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/genericPatchFields',
    '-I' + recursive_include_dirs / 'src/finiteArea',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/genericPatchFields',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_genericPatchFields = library(
    'genericPatchFields',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_genericPatchFields)
