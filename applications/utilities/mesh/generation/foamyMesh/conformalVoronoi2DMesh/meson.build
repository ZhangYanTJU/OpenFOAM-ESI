srcfiles = [
    lnInclude_hack,
    'cv2DControls/cv2DControls.C',
]
link_with = []
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'applications/utilities/mesh/generation/foamyMesh/conformalVoronoi2DMesh',
    '-I' + recursive_include_dirs / 'applications/utilities/mesh/generation/foamyMesh/conformalVoronoi2DMesh',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_cv2DMesh = library(
    'cv2DMesh',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_cv2DMesh)
