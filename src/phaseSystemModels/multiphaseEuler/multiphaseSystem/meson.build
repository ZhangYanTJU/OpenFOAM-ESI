srcfiles = [
    lnInclude_hack,
    'phaseModel/phaseModel.C',
    'multiphaseSystem/multiphaseSystem.C',
    'diameterModels/diameterModel/diameterModel.C',
    'diameterModels/constantDiameter/constantDiameter.C',
    'diameterModels/isothermalDiameter/isothermalDiameter.C',
    'interfacialModels/dragModels/dragModel/dragModel.C',
    'interfacialModels/dragModels/Ergun/Ergun.C',
    'interfacialModels/dragModels/GidaspowErgunWenYu/GidaspowErgunWenYu.C',
    'interfacialModels/dragModels/GidaspowSchillerNaumann/GidaspowSchillerNaumann.C',
    'interfacialModels/dragModels/SchillerNaumann/SchillerNaumann.C',
    'interfacialModels/dragModels/Gibilaro/Gibilaro.C',
    'interfacialModels/dragModels/WenYu/WenYu.C',
    'interfacialModels/dragModels/SyamlalOBrien/SyamlalOBrien.C',
    'interfacialModels/dragModels/blended/blended.C',
    'interfacialModels/dragModels/interface/interface.C',
    'interfacialModels/heatTransferModels/heatTransferModel/heatTransferModel.C',
    'interfacialModels/heatTransferModels/RanzMarshall/RanzMarshall.C',
    'derivedFvPatchFields/alphaContactAngle/alphaContactAngleFvPatchScalarField.C',
]
link_with = [
    lib_finiteVolume,
    lib_meshTools,
    lib_interfaceProperties,
    lib_incompressibleTransportModels,
]
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/phaseSystemModels/multiphaseEuler/multiphaseSystem',
    '-I' + meson.source_root() / 'src/transportModels',
    '-I' + recursive_include_dirs / 'src/transportModels/incompressible',
    '-I' + recursive_include_dirs / 'src/transportModels/interfaceProperties',
    '-I' + recursive_include_dirs / 'src/finiteVolume',
    '-I' + recursive_include_dirs / 'src/meshTools',
    '-I' + recursive_include_dirs / 'src/phaseSystemModels/multiphaseEuler/multiphaseSystem',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
lib_multiphaseSystem = library(
    'multiphaseSystem',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_multiphaseSystem)
