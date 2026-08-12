#!/bin/bash
# Check out a specific git commit in each module repository listed in a
# hashes file (as produced by parsing a bootstrap update log).
#
# Usage:
#   checkout_module_hashes.sh [-n|--dry-run] [hashes_file] [modules_dir]
#
# Defaults:
#   hashes_file: module_hashes.txt next to this script
#   modules_dir: parent directory of this script's directory
#
# File format: "<module_directory> <40-char hash>" per line; blank lines and
# text after '#' are ignored.  Repositories are left in detached-HEAD state.
# Repositories with local modifications are NOT forced; git will refuse and
# the module is reported as failed.

set -u

script_dir=$(cd "$(dirname "$0")" && pwd)
dry_run=0
args=()
for arg in "$@"; do
  case "$arg" in
    -n|--dry-run) dry_run=1 ;;
    -h|--help) sed -n '2,15p' "$0"; exit 0 ;;
    *) args+=("$arg") ;;
  esac
done
hashes_file=${args[0]:-$script_dir/module_hashes.txt}
modules_dir=${args[1]:-$(dirname "$script_dir")}

if [ ! -f "$hashes_file" ]; then
  echo "error: hashes file not found: $hashes_file" >&2
  exit 2
fi
if [ ! -d "$modules_dir" ]; then
  echo "error: modules directory not found: $modules_dir" >&2
  exit 2
fi

n_ok=0
n_skip=0
failed=()

while IFS= read -r line; do
  line=${line%%#*}                     # strip comments
  # shellcheck disable=SC2086
  set -- $line                         # split into fields
  [ $# -eq 0 ] && continue             # blank line
  if [ $# -ne 2 ]; then
    echo "WARN: malformed line ignored: $line" >&2
    continue
  fi
  module=$1
  hash=$2
  repo=$modules_dir/$module

  if [ ! -d "$repo" ]; then
    echo "SKIP  $module: directory not found"
    n_skip=$((n_skip + 1))
    continue
  fi
  if ! git -C "$repo" rev-parse --git-dir >/dev/null 2>&1; then
    echo "SKIP  $module: not a git repository"
    n_skip=$((n_skip + 1))
    continue
  fi

  current=$(git -C "$repo" rev-parse HEAD 2>/dev/null)
  if [ "$current" = "$hash" ]; then
    echo "OK    $module: already at $hash"
    n_ok=$((n_ok + 1))
    continue
  fi

  if [ "$dry_run" -eq 1 ]; then
    echo "DRY   $module: would check out $hash (currently ${current:-unknown})"
    continue
  fi

  # make sure the commit is present locally; fetch if it is not
  if ! git -C "$repo" rev-parse --verify --quiet "$hash^{commit}" >/dev/null; then
    echo "FETCH $module: $hash not present locally, fetching..."
    git -C "$repo" fetch --tags origin >/dev/null 2>&1
  fi
  if ! git -C "$repo" rev-parse --verify --quiet "$hash^{commit}" >/dev/null; then
    echo "FAIL  $module: commit $hash not found even after fetch"
    failed+=("$module")
    continue
  fi

  if git -C "$repo" checkout --detach "$hash" >/dev/null 2>&1; then
    echo "OK    $module: checked out $hash"
    n_ok=$((n_ok + 1))
  else
    echo "FAIL  $module: checkout refused (local changes?)"
    failed+=("$module")
  fi
done < "$hashes_file"

echo
echo "done: $n_ok ok, $n_skip skipped, ${#failed[@]} failed"
if [ ${#failed[@]} -gt 0 ]; then
  echo "failed modules: ${failed[*]}"
  exit 1
fi
