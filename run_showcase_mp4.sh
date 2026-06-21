#!/usr/bin/env bash
#
# Render every showcase scene to an MP4.
#
# For each scenes/showcase/*.json it:
#   1. renders the scene with PNG output enabled (16-bit per-channel PNG sequence),
#   2. encodes the PNG sequence into an H.265 MP4 via ./pngs2mp4.sh,
#   3. moves the result to mp4/<scene_name>.mp4.
#
# Usage:
#   ./run_showcase_mp4.sh                 # render all showcase scenes
#   FPS=12 ./run_showcase_mp4.sh          # set MP4 frame rate (default 12)
#   LL=1 ./run_showcase_mp4.sh            # lossless H.265 (large files)
#   KEEP_PNGS=1 ./run_showcase_mp4.sh     # keep the intermediate PNG frames
#   SCENES="scenes/showcase/showcase_cell8.json" ./run_showcase_mp4.sh   # subset
#   DEBUG=1 ./run_showcase_mp4.sh         # quick low-res preview render
#
# Environment knobs forwarded to the renderer:
#   DEBUG, SKIP_LOCKS, ALWAYS_BVH, NEVER_BVH, etc. (see cmd/photons4d/main.go)
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

BIN="./photons4d"
export FPS="${FPS:-12}"

mkdir -p gifs pngs raws mp4

# Build the renderer if it is missing or out of date.
if [ ! -x "$BIN" ]; then
  echo "[build] $BIN not found, building release binary..."
  go build -o "$BIN" ./cmd/photons4d
fi

if [ -z "${SCENES:-}" ]; then
  SCENES=$(find ./scenes/showcase -iname '*.json' | sort)
fi

if [ -z "$SCENES" ]; then
  echo "No showcase scenes found in scenes/showcase/*.json"
  exit 1
fi

export PNG=1

count=0
for f in $SCENES; do
  count=$((count + 1))
  scene_name=$(basename "${f%.json}")
  echo "================================================================"
  echo "[$count] Rendering: $f"
  echo "================================================================"

  # Extract gifOut from the scene to know the PNG prefix the renderer will use.
  gif_out=$(sed -n 's/^[[:space:]]*"gifOut"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "$f" | head -n 1)
  if [ -z "$gif_out" ]; then
    echo "Cannot extract gifOut from $f, skipping"
    continue
  fi
  prefix=$(basename "${gif_out%.gif}")
  # DEBUG mode appends _debug to the gif/png name (see internal/photons4d/run.go).
  if [ -n "${DEBUG:-}" ]; then
    prefix="${prefix}_debug"
  fi

  "$BIN" "$f"

  if compgen -G "pngs/${prefix}_*.png" >/dev/null; then
    ./pngs2mp4.sh "$prefix" </dev/null
    mv -f "${prefix}.mp4" "mp4/${scene_name}.mp4"
    echo "[ok] mp4/${scene_name}.mp4"
    if [ -z "${KEEP_PNGS:-}" ]; then
      rm -f "pngs/${prefix}"_*.png
    fi
  else
    echo "No PNGs produced for ${prefix}; skipping MP4 for ${scene_name}"
  fi
done

echo
echo "Done. Rendered $count scene(s). MP4s are in: $ROOT/mp4/"
ls -la mp4/*.mp4 2>/dev/null || true
