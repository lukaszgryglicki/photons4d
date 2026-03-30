These are GPT-retuned scene JSON variants derived from the uploaded photons4d scene set.

Naming:
- original: scene-name.json
- retuned:  scene-name-gpt.json

What was changed in general:
- gifOut renamed to ./gifs/<scene-name>-gpt.gif
- weak single-light scenes were given tighter key lights and one or two low-intensity fill/rim lights
- very low-SPP scenes were raised to reduce noise
- maxBounces and gamma were raised for scenes that looked too dim or too shallow
- escape / angular-bin scenes were normalized and, where useful, given denser angular resolutions
- malformed lowercase vector/color/light keys in experimental scenes were normalized to the schema the Go loader expects

A per-file summary is in MANIFEST.json.
