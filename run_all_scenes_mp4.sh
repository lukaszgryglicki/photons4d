#!/bin/bash
# PNG=1 is set here automatically.
# EXCLUDE_GLOB="ultimate_all_objects_trap.json" ./run_all_scenes_mp4.sh
# KEEP_PNGS=1 ./run_all_scenes_mp4.sh
set -e
export PNG=1
if [ -z "$SCENES" ]
then
  if [ -n "${EXCLUDE_GLOB:-}" ]; then
    SCENES=$(find ./scenes/ -iname "*.json" ! -iname "$EXCLUDE_GLOB")
  else
    SCENES=$(find ./scenes/ -iname "*.json")
  fi
fi
for f in $SCENES
do
  echo "Processing scene: $f"
  gif_out=$(sed -n 's/^[[:space:]]*"gifOut"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$f" | head -n 1)
  if [ -z "$gif_out" ]; then
    echo "Cannot extract gifOut from $f"
    exit 1
  fi
  prefix=$(basename "${gif_out%.gif}")
  scene_name=$(basename "${f%.json}")
  ./photons4d "$f"
  ./pngs2mp4.sh "$prefix" </dev/null
  if compgen -G "pngs/${prefix}_*.png" >/dev/null; then
    ./pngs2mp4.sh "$prefix" </dev/null
  else
    echo "No PNGs for ${prefix}, skipping MP4"
  fi
  mv -f "${prefix}.mp4" "mp4/${scene_name}.mp4"
  if [ -z "${KEEP_PNGS:-}" ]; then
    rm -f "pngs/${prefix}"_*.png
  fi
done
