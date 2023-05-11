srcfiles = [
    lnInclude_hack,
    'POSIX.C',
    'cpuInfo/cpuInfo.C',
    'cpuTime/cpuTimePosix.C',
    'memInfo/memInfo.C',
    'signals/sigFpe.C',
    'signals/sigSegv.C',
    'signals/sigInt.C',
    'signals/sigQuit.C',
    'signals/sigStopAtWriteNow.C',
    'signals/sigWriteNow.C',
    'signals/timer.C',
    'regExp/regExpPosix.C',
    'fileStat/fileStat.C',
    'fileMonitor/fileMonitor.C',
    'printStack/printStack.C',
]
link_with = []
dependencies = []
cpp_args = [
    '-I' + meson.source_root() / 'src/OSspecific/POSIX',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
    '-I' + recursive_include_dirs / 'src/OpenFOAM',
    '-I' + recursive_include_dirs / 'src/OSspecific/POSIX',
]
if fs.is_file('/usr/include/sys/inotify.h')
    cpp_args += '-DFOAM_USE_INOTIFY'
endif
lib_OSspecific = library(
    'OSspecific',
    srcfiles,
    link_with: link_with,
    dependencies: dependencies,
    cpp_args: cpp_args,
    implicit_include_directories: false,
    install: true,
)
pkg.generate(lib_OSspecific)
