#!/bin/bash
# Record the current git commit hash of every repository under a modules
# directory, in the format consumed by checkout_module_hashes.sh:
#   <module_directory> <commit_hash>  # <branch>[,dirty]
#
# Usage:
#   get_module_hashes.sh [modules_dir] > module_hashes.txt
#
# Default modules_dir: parent directory of this script's directory.
# Covers repositories one and two levels deep (e.g. chem_data/geostd).
# Submodules (.git is a file, state captured by the parent repo) are skipped.
# "dirty" marks repositories with uncommitted changes to tracked files.

set -u

script_dir=$(cd "$(dirname "$0")" && pwd)
modules_dir=${1:-$(dirname "$script_dir")}

if [ ! -d "$modules_dir" ]; then
  echo "error: modules directory not found: $modules_dir" >&2
  exit 2
fi
modules_dir=$(cd "$modules_dir" && pwd)

echo "# Module git hashes recorded from $modules_dir"
echo "# $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "# format: <module_directory> <commit_hash>  # <branch>[,dirty]"

find "$modules_dir" -mindepth 2 -maxdepth 3 -name .git -type d 2>/dev/null \
    | sort | while IFS= read -r gitpath; do
  repo=$(dirname "$gitpath")
  rel=${repo#"$modules_dir"/}
  if ! hash=$(git -C "$repo" rev-parse HEAD 2>/dev/null); then
    echo "WARN: $rel: cannot resolve HEAD (empty repository?), skipped" >&2
    continue
  fi
  branch=$(git -C "$repo" symbolic-ref --short -q HEAD)
  [ -z "$branch" ] && branch=detached
  dirty=""
  if [ -n "$(git -C "$repo" status --porcelain -uno 2>/dev/null)" ]; then
    dirty=",dirty"
  fi
  printf '%-20s %s  # %s%s\n' "$rel" "$hash" "$branch" "$dirty"
done
