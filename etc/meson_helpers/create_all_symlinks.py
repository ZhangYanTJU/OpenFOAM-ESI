#!/usr/bin/env python3
import sys
from os import path, listdir, walk
import os
from pathlib import Path

source_root = Path(sys.argv[1])
build_root = Path(sys.argv[2])


def create_symlinks_for_dir(subdir):
    outdir = build_root / subdir.relative_to(source_root)
    outdir.mkdir(parents=True, exist_ok=True)
    for fp in subdir.rglob("*.[CHh]"):
        if "lnInclude" in fp.parts:
            continue
        if (
            fp.name
            in [  # ugly name collisions. I hope this does not result in any problems.
                "fieldExprLemonParser.h",
                "patchExprLemonParser.h",
                "volumeExprLemonParser.h",
            ]
        ):
            continue
        outfile = outdir / fp.name
        if outfile.is_symlink():
            if os.readlink(outfile) != fp:
                outfile.unlink()
                outfile.symlink_to(fp)
        else:
            outfile.symlink_to(fp)


for subdir in source_root.rglob("Make"):
    if not path.isdir(subdir):
        continue
    subdir = subdir.parent
    create_symlinks_for_dir(subdir)

for (
    el
) in [  # These are the only directories found using `rg wmakeLnInclude` that do not have a `Make`` subdirectory
    "src/TurbulenceModels/phaseCompressible",
    "src/TurbulenceModels/phaseIncompressible",
]:
    create_symlinks_for_dir(source_root / el)

backlink = (build_root / "source")
if backlink.is_symlink():
    assert(os.readlink(backlink) == source_root)
else:
    backlink.symlink_to(source_root)

Path(build_root / "fake.h").touch()  # To make sure this script is not rerun nedlessly

# srcroot=$1
# recdir=$2
# targetdir=$PWD${recdir#"$srcroot"}

# for el in $(find "$recdir" -name '*.[CH]' ! -path '*lnInclude*'); do
#     # We use ${el##*/} instead of $(basename ${el}), because the latter is slower
#     target="$targetdir/${el##*/}"
#     if [ ! -e "$target" ]; then
#         ln -s "$el" "$target"
#     fi
# done

# # rootdir=$1
# # shift
# # for el in "$@"; do
# #     echo "$el"
# #     if [ ! -f "$el" ]; then
# #         ln -s "$rootdir/$el" "$el"
# #     fi
# # done
