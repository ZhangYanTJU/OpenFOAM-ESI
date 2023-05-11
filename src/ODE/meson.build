srcfiles = [
    lnInclude_hack,
    'ODESolvers/ODESolver/ODESolver.C',
    'ODESolvers/ODESolver/ODESolverNew.C',
    'ODESolvers/adaptiveSolver/adaptiveSolver.C',
    'ODESolvers/Euler/Euler.C',
    'ODESolvers/EulerSI/EulerSI.C',
    'ODESolvers/Trapezoid/Trapezoid.C',
    'ODESolvers/RKF45/RKF45.C',
    'ODESolvers/RKCK45/RKCK45.C',
    'ODESolvers/RKDP45/RKDP45.C',
    'ODESolvers/Rosenbrock12/Rosenbrock12.C',
    'ODESolvers/Rosenbrock23/Rosenbrock23.C',
    'ODESolvers/Rosenbrock34/Rosenbrock34.C',
    'ODESolvers/rodas23/rodas23.C',
    'ODESolvers/rodas34/rodas34.C',
    'ODESolvers/SIBS/SIBS.C',
    'ODESolvers/SIBS/SIMPR.C',
    'ODESolvers/SIBS/polyExtrapolate.C',
    'ODESolvers/seulex/seulex.C',
]
link_with = []
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/ODE',
    '-I' + recursive_include_dirs / 'src/ODE',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_ODE = library(
    'ODE',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_ODE)
