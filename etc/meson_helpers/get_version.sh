#!/usr/bin/env bash
root_path="$1"
source "$root_path/etc/bashrc"

set -euo pipefail
IFS=$'\n\t'
$root_path/wmake/scripts/wmake-build-info  | sed '1d;s/ = /=/g;s/    //g;s/\n/kk/g' | tr '\n' ';' | sed 's/;/; /g;s/; $/\n/'
