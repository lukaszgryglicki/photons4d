#!/usr/bin/env bash
# THREADS=512 ALWAYS_BVH=1 ./run_all_scenes_mp4.sh
# FPS=5 ./run_all_scenes_mp4.sh
# LOSSLESS=1 ./run_all_scenes_mp4.sh
# EXCLUDE_GLOB="ultimate_all_objects_trap.json" ./run_all_scenes_mp4.sh
set -euo pipefail

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
BIN="${BIN:-${ROOT}/photons4d}"
FPS="${FPS:-}"                 # optional override; if empty, derive from gifDelay
LOSSLESS="${LOSSLESS:-0}"      # 1 => x265 lossless
THREADS="${THREADS:-}"         # optional; exported to GOMAXPROCS if set
EXCLUDE_GLOB="${EXCLUDE_GLOB:-}" # e.g. ultimate_all_objects_trap.json

mkdir -p "${ROOT}/gifs" "${ROOT}/pngs" "${ROOT}/mp4" "${ROOT}/logs"

if [[ ! -x "${BIN}" ]]; then
  echo "Building ${BIN}"
  (cd "${ROOT}" && go build -trimpath -o "${BIN}" ./cmd/photons4d)
fi

if [[ -n "${THREADS}" ]]; then
  export GOMAXPROCS="${THREADS}"
fi

# Required for PNG frame generation.
export PNG=1

# Leave these untouched so the caller can still do:
#   ALWAYS_BVH=1 FORCE_ESCAPE=1 SPP_ADJUST=1 RAW=1 ./run_all_scenes_mp4.sh

status=0

while IFS= read -r -d '' scene; do
  base="$(basename "${scene}")"
  if [[ -n "${EXCLUDE_GLOB}" && "${base}" == ${EXCLUDE_GLOB} ]]; then
    echo "Skipping scene: ${scene}"
    continue
  fi

  scene_id="${base%.json}"
  echo "Processing scene: ${scene}"

  mapfile -t meta < <(python3 - "$scene" <<'PY'
import json, sys

path = sys.argv[1]
with open(path, "r", encoding="utf-8") as fh:
    cfg = json.load(fh)

gif_out = cfg.get("gifOut") or "./gifs/volume.gif"
gif_delay = cfg.get("gifDelay") or 5

if gif_out.lower().endswith(".gif"):
    png_prefix = gif_out[:-4]
else:
    png_prefix = gif_out

png_prefix = png_prefix.replace("gifs/", "pngs/", 1)

fps = 100.0 / float(gif_delay)

print(gif_out)
print(png_prefix)
print(f"{fps:.6f}".rstrip("0").rstrip("."))
PY
  )

  gif_out="${meta[0]}"
  png_prefix="${meta[1]}"
  derived_fps="${meta[2]}"

  fps_use="${FPS:-${derived_fps}}"
  mp4_out="${ROOT}/mp4/${scene_id}.mp4"
  log_out="${ROOT}/logs/${scene_id}_$(date +%Y%m%d_%H%M%S).log"

  mkdir -p "$(dirname "${mp4_out}")" "$(dirname "${png_prefix}")"

  # Remove stale frames from an older failed/aborted run of the same output prefix.
  rm -f "${png_prefix}"_*.png

  if ! "${BIN}" "${scene}" 2>&1 | tee "${log_out}"; then
    echo "Render failed for ${scene}; leaving any partial PNGs for inspection." >&2
    status=1
    continue
  fi

  if ! compgen -G "${png_prefix}_*.png" >/dev/null; then
    echo "No PNG frames found for ${scene} under ${png_prefix}_*.png" >&2
    status=1
    continue
  fi

  if [[ "${LOSSLESS}" == "1" ]]; then
    if ! ffmpeg -y -framerate "${fps_use}" -pattern_type glob -i "${png_prefix}_*.png" \
      -c:v libx265 -preset slow -pix_fmt yuv444p12le \
      -x265-params lossless=1 \
      -movflags +faststart "${mp4_out}"; then
      echo "ffmpeg failed for ${scene}; keeping PNGs under ${png_prefix}_*.png" >&2
      status=1
      continue
    fi
  else
    if ! ffmpeg -y -framerate "${fps_use}" -pattern_type glob -i "${png_prefix}_*.png" \
      -c:v libx265 -preset slow -pix_fmt yuv444p12le \
      -x265-params crf=2:aq-mode=3:aq-strength=1.2:psy-rd=2.0:psy-rdoq=1.5 \
      -tune grain \
      -movflags +faststart "${mp4_out}"; then
      echo "ffmpeg failed for ${scene}; keeping PNGs under ${png_prefix}_*.png" >&2
      status=1
      continue
    fi
  fi

  rm -f "${png_prefix}"_*.png

  echo "GIF : ${gif_out}"
  echo "MP4 : ${mp4_out}"
  echo "LOG : ${log_out}"
  echo
done < <(find "${ROOT}/scenes" -type f -iname "*.json" -print0)

exit "${status}"
