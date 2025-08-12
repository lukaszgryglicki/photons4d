#!/bin/bash
# LL=1 # Set to 1 for lossless encoding
if [ -z "$1" ]
then
  echo "Usage: $0 <prefix>"
  exit 1
fi
prefix="$1"
if [ -z "$LL" ]
then
  ffmpeg -y -framerate 24 -pattern_type glob -i "pngs/${prefix}_*.png" -c:v libx265 -preset slow -pix_fmt yuv444p12le -x265-params crf=2:tune=grain:aq-mode=3 -movflags +faststart "${prefix}.mp4"
else
  ffmpeg -y -framerate 24 -pattern_type glob -i "pngs/${prefix}_*.png" -c:v libx265 -preset slow -pix_fmt yuv444p12le -x265-params lossless=1 -movflags +faststart "${prefix}.mp4"
fi
