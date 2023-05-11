srcfiles = [
    lnInclude_hack,
    'equilibriumFlameT.C',
]
link_with = [
    lib_OpenFOAM,
    lib_specie,
]
dependencies = [
    m_dep,
    dl_dep,
]
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/thermophysical/equilibriumFlameT',
    '-I' + recursive_include_dirs / 'src/thermophysicalModels/specie',
    '-I' + recursive_include_dirs / 'applications/utilities/thermophysical/equilibriumFlameT',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
exe_equilibriumFlameT = executable(
    'equilibriumFlameT',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
