srcfiles = [
    lnInclude_hack,
    'potential/potential.C',
    'pairPotential/pairPotentialList/pairPotentialList.C',
    'pairPotential/basic/pairPotential.C',
    'pairPotential/basic/pairPotentialIO.C',
    'pairPotential/basic/pairPotentialNew.C',
    'pairPotential/derived/lennardJones/lennardJones.C',
    'pairPotential/derived/maitlandSmith/maitlandSmith.C',
    'pairPotential/derived/azizChen/azizChen.C',
    'pairPotential/derived/exponentialRepulsion/exponentialRepulsion.C',
    'pairPotential/derived/coulomb/coulomb.C',
    'pairPotential/derived/dampedCoulomb/dampedCoulomb.C',
    'pairPotential/derived/noInteraction/noInteraction.C',
    'energyScalingFunction/basic/energyScalingFunction.C',
    'energyScalingFunction/basic/energyScalingFunctionNew.C',
    'energyScalingFunction/derived/shifted/shifted.C',
    'energyScalingFunction/derived/shiftedForce/shiftedForce.C',
    'energyScalingFunction/derived/noScaling/noScaling.C',
    'energyScalingFunction/derived/sigmoid/sigmoid.C',
    'energyScalingFunction/derived/doubleSigmoid/doubleSigmoid.C',
    'tetherPotential/tetherPotentialList/tetherPotentialList.C',
    'tetherPotential/basic/tetherPotential.C',
    'tetherPotential/basic/tetherPotentialNew.C',
    'tetherPotential/derived/harmonicSpring/harmonicSpring.C',
    'tetherPotential/derived/restrainedHarmonicSpring/restrainedHarmonicSpring.C',
    'tetherPotential/derived/pitchForkRing/pitchForkRing.C',
    'electrostaticPotential/electrostaticPotential.C',
]
link_with = [
    lib_finiteVolume,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/lagrangian/molecularDynamics/potential',
    '-I' + meson.source_root() / 'src/lagrangian/molecularDynamics',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/lagrangian/molecularDynamics/potential',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_potential = library(
    'potential',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_potential)
