#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-$(pwd)}"
SCENE_NAME="ultimate_all_objects_trap"
CFG="${ROOT}/scenes/${SCENE_NAME}.json"
BIN="${ROOT}/photons4d"
THREADS="${THREADS:-512}"
FPS="${FPS:-5}"
LOSSLESS="${LOSSLESS:-0}"

mkdir -p "${ROOT}/scenes" "${ROOT}/gifs" "${ROOT}/pngs" "${ROOT}/raws" "${ROOT}/mp4" "${ROOT}/logs"

export CFG
python3 - <<'PY'
import json, math, os

def p4(x,y,z,w): return {"X": round(x,4), "Y": round(y,4), "Z": round(z,4), "W": round(w,4)}
def v4(x,y,z,w): return {"X": round(x,4), "Y": round(y,4), "Z": round(z,4), "W": round(w,4)}
def rot(xy,xz,xw,yz,yw,zw): return {"xy": xy, "xz": xz, "xw": xw, "yz": yz, "yw": yw, "zw": zw}
def rgb(r,g,b): return {"R": round(r,4), "G": round(g,4), "B": round(b,4)}
def mat(color, diffuse, reflect, refract, ior):
    return {
        "color": rgb(*color),
        "diffuse": rgb(*diffuse),
        "reflect": rgb(*reflect),
        "refract": rgb(*refract),
        "ior": rgb(*ior),
    }

cfg = {
    "sceneResX": 256,
    "sceneResY": 256,
    "sceneResZ": 16,
    "probeRays": 65536,
    "spp": 256,
    "gifOut": "./gifs/ultimate_all_objects_trap.gif",
    "gifDelay": 20,
    "gamma": 0.8,
    "scene": {
        "center": {"X": 0, "Y": 0, "Z": 0, "W": 0},
        "width": 2.2,
        "height": 2.2,
        "depth": 2.2,
        "maxBounces": 64
    },
    "lights": [],
    "hyperspheres": [],
    "cells5": [],
    "cells8": [],
    "cells16": [],
    "cells24": [],
    "cells120": [],
    "cells600": [],
}

cfg["lights"] = [
    {"origin": p4(0.0, 0.0, 0.0, 1.16), "direction": v4(0.0, 0.0, 0.0, -1.0), "color": rgb(1.0, 1.0, 1.0), "angleDeg": 44.0, "intensity": 0.95},
    {"origin": p4(0.92, 0.0, 0.0, 1.08), "direction": v4(-0.42, 0.0, 0.0, -1.0), "color": rgb(1.0, 0.35, 0.30), "angleDeg": 24.0, "intensity": 0.60},
    {"origin": p4(-0.92, 0.0, 0.0, 1.08), "direction": v4(0.42, 0.0, 0.0, -1.0), "color": rgb(0.35, 0.95, 1.0), "angleDeg": 24.0, "intensity": 0.60},
    {"origin": p4(0.0, 0.92, 0.0, 1.08), "direction": v4(0.0, -0.42, 0.0, -1.0), "color": rgb(0.40, 1.0, 0.45), "angleDeg": 24.0, "intensity": 0.58},
    {"origin": p4(0.0, -0.92, 0.0, 1.08), "direction": v4(0.0, 0.42, 0.0, -1.0), "color": rgb(1.0, 0.90, 0.35), "angleDeg": 24.0, "intensity": 0.58},
    {"origin": p4(0.72, 0.72, 0.55, 1.12), "direction": v4(-0.35, -0.35, -0.20, -1.0), "color": rgb(0.95, 0.45, 1.0), "angleDeg": 18.0, "intensity": 0.44},
    {"origin": p4(-0.72, 0.72, -0.55, 1.12), "direction": v4(0.35, -0.35, 0.20, -1.0), "color": rgb(0.30, 0.75, 1.0), "angleDeg": 18.0, "intensity": 0.44},
    {"origin": p4(0.72, -0.72, -0.55, 1.12), "direction": v4(-0.35, 0.35, 0.20, -1.0), "color": rgb(1.0, 0.55, 0.35), "angleDeg": 18.0, "intensity": 0.44},
    {"origin": p4(-0.72, -0.72, 0.55, 1.12), "direction": v4(0.35, 0.35, -0.20, -1.0), "color": rgb(0.55, 1.0, 0.70), "angleDeg": 18.0, "intensity": 0.44},
    {"origin": p4(0.22, -0.18, 0.10, 0.86), "direction": v4(-0.20, 0.12, -0.05, -1.0), "color": rgb(1.0, 0.55, 0.20), "angleDeg": 10.0, "intensity": 0.24},
    {"origin": p4(-0.26, 0.22, -0.12, 0.84), "direction": v4(0.18, -0.10, 0.07, -1.0), "color": rgb(0.25, 0.55, 1.0), "angleDeg": 10.0, "intensity": 0.24},
    {"origin": p4(0.0, -0.12, 0.0, 0.95), "direction": v4(0.0, 0.08, 0.0, -1.0), "color": rgb(0.22, 0.22, 0.28), "angleDeg": 30.0, "intensity": -0.14},
]

cfg["cells8"].append({
    "center": p4(0.0, 0.0, 0.0, 0.56),
    "scale": v4(1.72, 1.72, 1.72, 0.76),
    "rotDeg": rot(12, -8, 9, 6, -10, 14),
    **mat((0.96, 0.98, 1.0), (0.03, 0.03, 0.03), (0.58, 0.62, 0.66), (0.34, 0.30, 0.26), (1.42, 1.46, 1.50))
})
cfg["hyperspheres"].append({
    "center": p4(0.02, -0.03, 0.04, 0.58),
    "scale": v4(0.82, 0.72, 0.76, 0.24),
    "rotDeg": rot(22, 14, -10, 26, 8, -18),
    **mat((0.96, 0.98, 1.0), (0.02, 0.02, 0.02), (0.08, 0.08, 0.08), (0.88, 0.86, 0.84), (1.60, 1.48, 1.72))
})
cfg["cells120"].append({
    "center": p4(0.0, 0.0, 0.0, 0.60),
    "scale": v4(0.72, 0.72, 0.72, 0.72),
    "rotDeg": rot(28, -16, 12, 10, -18, 24),
    **mat((0.90, 0.95, 1.0), (0.04, 0.04, 0.04), (0.30, 0.34, 0.38), (0.64, 0.60, 0.56), (1.46, 1.38, 1.56))
})
cfg["cells600"].append({
    "center": p4(0.0, 0.0, 0.0, 0.62),
    "scale": v4(0.58, 0.58, 0.58, 0.58),
    "rotDeg": rot(-18, 24, -12, 16, 22, -14),
    **mat((1.0, 0.96, 0.92), (0.08, 0.08, 0.08), (0.70, 0.64, 0.58), (0.18, 0.24, 0.30), (1.18, 1.24, 1.30))
})

angles = [15, 78, 141, 204, 267, 330]

mats5 = [
    ((1.0,0.35,0.30),(0.72,0.62,0.62),(0.15,0.14,0.14),(0.08,0.14,0.18),(1.08,1.10,1.14)),
    ((0.90,0.95,1.0),(0.08,0.08,0.08),(0.22,0.22,0.22),(0.66,0.64,0.62),(1.56,1.40,1.68)),
    ((0.35,0.85,1.0),(0.42,0.42,0.42),(0.32,0.32,0.32),(0.21,0.21,0.21),(1.26,1.34,1.44)),
    ((1.0,0.92,0.45),(0.22,0.16,0.12),(0.56,0.52,0.44),(0.16,0.22,0.34),(1.34,1.50,1.72)),
    ((0.42,1.0,0.45),(0.62,0.62,0.62),(0.20,0.20,0.20),(0.12,0.12,0.12),(1.12,1.12,1.12)),
    ((0.95,0.45,1.0),(0.10,0.10,0.10),(0.18,0.26,0.34),(0.67,0.57,0.47),(1.46,1.34,1.62)),
]
for idx,a in enumerate(angles):
    th = math.radians(a)
    r = 0.63
    x = r*math.cos(th)
    y = r*math.sin(th)
    z = 0.45*math.sin(2*th)
    w = 0.69 + 0.06*math.cos(1.5*th)
    scale = (0.10 + 0.02*(idx%2), 0.08 + 0.03*((idx+1)%3==0), 0.08 + 0.02*(idx%3==1), 0.08 + 0.03*(idx%2))
    color,diff,ref,refr,ior = mats5[idx]
    cfg["cells5"].append({"center": p4(x,y,z,w), "scale": v4(*scale), "rotDeg": rot((a*1.3)%60-30,(a*0.9)%70-35,(a*1.7)%80-40,(a*1.1)%60-30,(a*0.7)%70-35,(a*1.5)%80-40), **mat(color,diff,ref,refr,ior)})

mats8 = [
    ((1.0,0.25,0.25),(0.82,0.82,0.82),(0.10,0.10,0.10),(0.03,0.03,0.03),(1.0,1.0,1.0)),
    ((0.25,0.65,1.0),(0.12,0.12,0.12),(0.78,0.78,0.78),(0.05,0.05,0.05),(1.0,1.0,1.0)),
    ((0.95,0.95,1.0),(0.03,0.03,0.03),(0.10,0.10,0.10),(0.86,0.84,0.82),(1.62,1.46,1.56)),
    ((0.30,1.0,0.35),(0.56,0.56,0.56),(0.18,0.12,0.18),(0.20,0.26,0.20),(1.34,1.54,1.24)),
    ((1.0,0.85,0.25),(0.0,0.0,0.0),(0.92,0.82,0.62),(0.04,0.08,0.18),(1.0,1.0,1.0)),
    ((0.22,0.90,1.0),(0.20,0.20,0.20),(0.34,0.34,0.34),(0.42,0.42,0.42),(1.46,1.52,1.38)),
]
for idx,a in enumerate([0,60,120,180,240,300]):
    th = math.radians(a+18)
    r = 0.44
    x = r*math.cos(th)
    y = r*math.sin(th)
    z = 0.32*math.cos(2*th)
    w = 0.56 + 0.05*math.sin(1.2*th)
    scale = (0.16+0.02*(idx%3==0), 0.14+0.04*(idx%2), 0.14+0.03*(idx%3==1), 0.12+0.02*(idx%2))
    color,diff,ref,refr,ior = mats8[idx]
    cfg["cells8"].append({"center": p4(x,y,z,w), "scale": v4(*scale), "rotDeg": rot((a*1.7)%80-40,(a*0.6)%70-35,(a*1.1)%60-30,(a*0.8)%80-40,(a*1.4)%70-35,(a*1.2)%80-40), **mat(color,diff,ref,refr,ior)})

mats16 = [
    ((0.95,0.65,0.65),(0.74,0.74,0.74),(0.12,0.08,0.12),(0.06,0.10,0.06),(1.0,1.0,1.0)),
    ((0.85,0.92,1.0),(0.08,0.08,0.08),(0.20,0.20,0.20),(0.70,0.70,0.70),(1.58,1.36,1.72)),
    ((0.25,0.25,1.0),(0.46,0.46,0.46),(0.38,0.38,0.38),(0.10,0.10,0.10),(1.24,1.24,1.24)),
    ((0.25,1.0,0.25),(0.0,0.0,0.0),(0.0,0.0,0.0),(0.96,0.96,0.96),(1.74,1.22,1.48)),
    ((1.0,0.42,0.15),(0.28,0.22,0.16),(0.58,0.52,0.44),(0.10,0.16,0.26),(1.0,1.0,1.0)),
    ((0.92,0.92,0.92),(0.18,0.18,0.18),(0.26,0.26,0.26),(0.50,0.50,0.50),(1.48,1.38,1.58)),
]
for idx,a in enumerate([30,90,150,210,270,330]):
    th = math.radians(a)
    r = 0.52
    x = r*math.cos(th)
    y = r*math.sin(th)
    z = 0.50*math.sin(1.5*th)
    w = 0.73 + 0.05*math.cos(2.2*th)
    scale = (0.12+0.02*(idx%2), 0.10+0.03*(idx%3==0), 0.10+0.03*(idx%3==1), 0.10+0.02*(idx%2))
    color,diff,ref,refr,ior = mats16[idx]
    cfg["cells16"].append({"center": p4(x,y,z,w), "scale": v4(*scale), "rotDeg": rot((a*0.8)%80-40,(a*1.2)%90-45,(a*0.5)%70-35,(a*1.1)%80-40,(a*0.9)%90-45,(a*1.4)%80-40), **mat(color,diff,ref,refr,ior)})

mats24 = [
    ((0.45,1.0,0.90),(0.58,0.58,0.58),(0.12,0.18,0.22),(0.16,0.16,0.10),(1.28,1.54,1.22)),
    ((1.0,0.32,0.62),(0.08,0.08,0.08),(0.82,0.72,0.62),(0.05,0.05,0.05),(1.0,1.0,1.0)),
    ((0.92,0.92,1.0),(0.05,0.05,0.05),(0.14,0.14,0.14),(0.76,0.76,0.76),(1.62,1.42,1.48)),
    ((0.22,0.62,1.0),(0.30,0.22,0.22),(0.50,0.42,0.42),(0.15,0.25,0.25),(1.28,1.38,1.34)),
    ((1.0,1.0,0.32),(0.62,0.62,0.12),(0.20,0.20,0.10),(0.12,0.12,0.12),(1.10,1.10,1.10)),
    ((0.92,0.92,0.92),(0.14,0.14,0.14),(0.16,0.20,0.16),(0.66,0.62,0.66),(1.52,1.42,1.58)),
]
for idx,a in enumerate([45,105,165,225,285,345]):
    th = math.radians(a)
    r = 0.36
    x = r*math.cos(th)
    y = r*math.sin(th)
    z = 0.54*math.cos(1.7*th)
    w = 0.63 + 0.06*math.sin(1.9*th)
    scale = (0.10+0.03*(idx%2), 0.08+0.03*(idx%3==1), 0.10+0.02*(idx%3==2), 0.08+0.03*(idx%2))
    color,diff,ref,refr,ior = mats24[idx]
    cfg["cells24"].append({"center": p4(x,y,z,w), "scale": v4(*scale), "rotDeg": rot((a*1.0)%70-35,(a*1.6)%90-45,(a*1.3)%80-40,(a*0.9)%70-35,(a*0.7)%90-45,(a*1.1)%80-40), **mat(color,diff,ref,refr,ior)})

matsHS = [
    ((1.0,1.0,1.0),(0.0,0.0,0.0),(0.92,0.92,0.92),(0.03,0.03,0.03),(1.0,1.0,1.0)),
    ((0.92,0.96,1.0),(0.08,0.08,0.08),(0.12,0.12,0.12),(0.76,0.76,0.76),(1.54,1.42,1.66)),
    ((1.0,0.35,0.35),(0.54,0.18,0.18),(0.24,0.20,0.20),(0.16,0.54,0.54),(1.28,1.58,1.58)),
    ((0.35,1.0,0.9),(0.16,0.26,0.26),(0.42,0.26,0.26),(0.36,0.42,0.42),(1.32,1.42,1.28)),
    ((1.0,0.92,0.55),(0.44,0.40,0.24),(0.28,0.28,0.24),(0.18,0.22,0.46),(1.22,1.38,1.62)),
]
for idx,a in enumerate([20,92,164,236,308]):
    th = math.radians(a)
    r = 0.72
    x = r*math.cos(th)
    y = r*math.sin(th)
    z = 0.38*math.sin(1.4*th)
    w = 0.78 + 0.04*math.cos(2.5*th)
    scale = (0.12-0.02*(idx%2), 0.10-0.02*(idx%3==1), 0.10+0.03*(idx%2), 0.08+0.02*(idx%3==2))
    color,diff,ref,refr,ior = matsHS[idx]
    cfg["hyperspheres"].append({"center": p4(x,y,z,w), "scale": v4(*scale), "rotDeg": rot((a*1.1)%90-45,(a*0.8)%70-35,(a*1.3)%90-45,(a*0.9)%70-35,(a*1.5)%90-45,(a*1.2)%70-35), **mat(color,diff,ref,refr,ior)})

mats120 = [
    ((1.0,0.35,0.35),(0.62,0.54,0.54),(0.18,0.10,0.10),(0.12,0.18,0.24),(1.34,1.44,1.56)),
    ((0.25,0.70,1.0),(0.16,0.16,0.16),(0.72,0.72,0.72),(0.07,0.07,0.07),(1.0,1.0,1.0)),
    ((0.96,0.96,1.0),(0.05,0.05,0.05),(0.12,0.12,0.12),(0.78,0.76,0.78),(1.60,1.42,1.52)),
    ((0.30,1.0,0.30),(0.54,0.54,0.54),(0.18,0.10,0.18),(0.18,0.26,0.18),(1.34,1.54,1.22)),
    ((1.0,0.84,0.25),(0.0,0.0,0.0),(0.92,0.82,0.58),(0.04,0.08,0.20),(1.0,1.0,1.0)),
]
for idx,a in enumerate([12,84,156,228,300]):
    th = math.radians(a)
    r = 0.56
    x = r*math.cos(th)
    y = r*math.sin(th)
    z = 0.28*math.cos(2*th)
    w = 0.60 + 0.04*math.sin(2.3*th)
    scale = (0.12+0.01*(idx%2), 0.12+0.01*(idx%3==0), 0.12+0.01*(idx%3==1), 0.12+0.01*(idx%2))
    color,diff,ref,refr,ior = mats120[idx]
    cfg["cells120"].append({"center": p4(x,y,z,w), "scale": v4(*scale), "rotDeg": rot((a*1.4)%90-45,(a*0.9)%70-35,(a*1.1)%80-40,(a*0.8)%90-45,(a*1.5)%70-35,(a*1.2)%80-40), **mat(color,diff,ref,refr,ior)})

mats600 = [
    ((0.45,1.0,0.90),(0.56,0.56,0.56),(0.12,0.18,0.22),(0.16,0.16,0.10),(1.28,1.54,1.22)),
    ((1.0,0.32,0.62),(0.08,0.08,0.08),(0.78,0.68,0.58),(0.10,0.18,0.18),(1.0,1.0,1.0)),
    ((0.92,0.92,1.0),(0.05,0.05,0.05),(0.14,0.14,0.14),(0.76,0.76,0.76),(1.62,1.42,1.48)),
    ((0.22,0.62,1.0),(0.28,0.22,0.22),(0.52,0.44,0.44),(0.15,0.24,0.24),(1.28,1.38,1.34)),
    ((1.0,1.0,0.32),(0.62,0.62,0.12),(0.20,0.20,0.10),(0.12,0.12,0.12),(1.10,1.10,1.10)),
]
for idx,a in enumerate([48,120,192,264,336]):
    th = math.radians(a)
    r = 0.48
    x = r*math.cos(th)
    y = r*math.sin(th)
    z = 0.26*math.sin(2.1*th)
    w = 0.66 + 0.05*math.cos(1.7*th)
    scale = (0.11+0.01*(idx%2), 0.11+0.01*(idx%3==1), 0.11+0.01*(idx%3==2), 0.11+0.01*(idx%2))
    color,diff,ref,refr,ior = mats600[idx]
    cfg["cells600"].append({"center": p4(x,y,z,w), "scale": v4(*scale), "rotDeg": rot((a*1.0)%90-45,(a*1.3)%70-35,(a*0.7)%80-40,(a*1.5)%90-45,(a*0.9)%70-35,(a*1.1)%80-40), **mat(color,diff,ref,refr,ior)})

out = os.environ["CFG"]
with open(out, "w", encoding="utf-8") as fh:
    json.dump(cfg, fh, indent=2)
    fh.write("\n")

print(out)
print("lights=", len(cfg["lights"]))
print("objects=", sum(len(cfg[k]) for k in ("hyperspheres","cells5","cells8","cells16","cells24","cells120","cells600")))
PY

go build -trimpath -o "${BIN}" ./cmd/photons4d

export GOMAXPROCS="${THREADS}"
export ALWAYS_BVH=1
export PNG=1
export RAW=1
unset SKIP_LOCKS
unset FORCE_ESCAPE
unset SPP_ADJUST
unset DEBUG

LOG="${ROOT}/logs/${SCENE_NAME}_$(date +%Y%m%d_%H%M%S).log"

time "${BIN}" "${CFG}" 2>&1 | tee "${LOG}"

if compgen -G "${ROOT}/pngs/${SCENE_NAME}_*.png" >/dev/null; then
  if [[ "${LOSSLESS}" == "1" ]]; then
    ffmpeg -y -framerate "${FPS}" -pattern_type glob -i "${ROOT}/pngs/${SCENE_NAME}_*.png" \
      -c:v libx265 -preset slow -pix_fmt yuv444p12le \
      -x265-params lossless=1 \
      -movflags +faststart "${ROOT}/mp4/${SCENE_NAME}.mp4"
  else
    ffmpeg -y -framerate "${FPS}" -pattern_type glob -i "${ROOT}/pngs/${SCENE_NAME}_*.png" \
      -c:v libx265 -preset slow -pix_fmt yuv444p12le \
      -x265-params crf=2:aq-mode=3:aq-strength=1.2:psy-rd=2.0:psy-rdoq=1.5 -tune grain \
      -movflags +faststart "${ROOT}/mp4/${SCENE_NAME}.mp4"
  fi
else
  echo "No PNG frames found under ${ROOT}/pngs/${SCENE_NAME}_*.png" >&2
  exit 1
fi

echo
echo "Done."
echo "GIF : ${ROOT}/gifs/${SCENE_NAME}.gif"
echo "RAW : ${ROOT}/raws/${SCENE_NAME}.raw"
echo "MP4 : ${ROOT}/mp4/${SCENE_NAME}.mp4"
echo "LOG : ${LOG}"
