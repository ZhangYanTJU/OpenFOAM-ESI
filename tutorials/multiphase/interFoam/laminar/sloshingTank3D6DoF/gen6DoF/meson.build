srcfiles = [
    lnInclude_hack,
    'gen6DoF.C',
]
link_with = [
    lib_OpenFOAM,
    lib_finiteVolume,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'tutorials/multiphase/interFoam/laminar/sloshingTank3D6DoF/gen6DoF',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'tutorials/multiphase/interFoam/laminar/sloshingTank3D6DoF/gen6DoF',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_gen6DoF = executable(
    'gen6DoF',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
    build_by_default: false,
)
