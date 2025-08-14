#!/bin/bash
# FORCE_ESCAPE=1 SPP_ADJUST=1 ./run_all_scenes.sh
rm gifs/*.gif
for f in $(find ./scenes/ -iname "*.json")
do
  echo "Processing scene: $f"
  ./photons4d "${f}"
done
