#!/usr/bin/env bash
root_path="$1"
output="$2"
shift
shift
source "$root_path/etc/bashrc"

set -euo pipefail
IFS=$'\n\t'

"$root_path/wmake/scripts/wmake-build-info" -update -filter "$root_path/src/OpenFOAM/global/foamConfig.Cver" > "$output"
