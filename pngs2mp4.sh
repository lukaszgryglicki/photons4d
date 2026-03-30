#!/bin/bash
# LL=1 FPS=5 ./pngs2mp4.sh prefix
if [ -z "$1" ]
then
  echo "Usage: $0 <prefix>"
  echo "Example: $0 120-cell-flat"
  exit 1
fi
prefix="$1"
fps="${FPS:-10}"
if [ -z "$LL" ]
then
  ffmpeg -nostdin -y -framerate "$fps" -pattern_type glob -i "pngs/${prefix}_*.png" -c:v libx265 -preset slow -pix_fmt yuv444p12le -x265-params crf=2:aq-mode=3:aq-strength=1.2:psy-rd=2.0:psy-rdoq=1.5 -tune grain -movflags +faststart "${prefix}.mp4"
else
  ffmpeg -nostdin -y -framerate "$fps" -pattern_type glob -i "pngs/${prefix}_*.png" -c:v libx265 -preset slow -pix_fmt yuv444p12le -x265-params lossless=1 -movflags +faststart "${prefix}.mp4"
fi
