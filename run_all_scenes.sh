#!/bin/bash
rm gifs/*.gif
for f in $(find ./scenes/ -iname "*.json")
do
  echo "Processing scene: $f"
  ./photons4d "${f}"
done
