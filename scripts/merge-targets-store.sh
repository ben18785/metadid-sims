#!/usr/bin/env bash
#
# Merge the per-scenario _targets store slices that the simulate jobs upload
# into the single store the report reads.
#
#   scripts/merge-targets-store.sh <artifacts-dir> [store-dir]
#
# Shared by sims.yaml and render-report.yaml so the two cannot drift. Target
# names are unique per category, so cross-category collisions only happen on
# shared dependencies like `metadid_src`, which all carry the same hash.
set -euo pipefail

artifacts=${1:?usage: merge-targets-store.sh <artifacts-dir> [store-dir]}
store=${2:-_targets}

shopt -s nullglob
slices=("$artifacts"/targets-store-*)
if [ ${#slices[@]} -eq 0 ]; then
  echo "::error::No targets-store-* slices found in $artifacts"
  exit 1
fi

mkdir -p "$store/meta" "$store/objects"
for d in "${slices[@]}"; do
  if [ -d "$d/objects" ]; then
    cp -rn "$d"/objects/. "$store/objects/" || true
  fi
done

metas=("$artifacts"/targets-store-*/meta/meta)
if [ ${#metas[@]} -eq 0 ]; then
  echo "::error::No meta files found in $artifacts"
  exit 1
fi

# Keep one header line, then dedupe rows by the first column (target name).
head -n 1 "${metas[0]}" > "$store/meta/meta"
for f in "${metas[@]}"; do
  tail -n +2 "$f"
done | awk -F'|' '!seen[$1]++' >> "$store/meta/meta"

echo "Merged store contains $(wc -l < "$store/meta/meta") meta rows and $(ls "$store/objects" | wc -l) object files."
