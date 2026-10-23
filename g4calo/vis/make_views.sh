#!/bin/bash
# 3D views for the project page: one event per sample drawn by Geant4 (VRML2FILE), converted to site/g4vis/*.json.
# The events are not selected: the first event of each run with a fixed seed.
set -e
cd "$(dirname "$0")/.."
source "$HOME/opt/geant4/bin/geant4.sh"
OUT=../site/g4vis
TMP=$(mktemp -d)
mkdir -p "$OUT"
while read -r p E X label; do
  name="${p}_${E}GeV_X0${X}"
  rm -f "$TMP"/g4_*.wrl
  G4VRMLFILE_DEST_DIR="$TMP/" build/calo "$p" "$E" 1 "$X" "$TMP/$name.txt" 7 macros/vrml.mac > /dev/null 2>&1
  python3 vis/wrl2json.py "$TMP/g4_00.wrl" "$OUT/$name.json" --label "$label"
done <<LIST
e- 50 2.0 e⁻, 50 GeV, 2 X₀ upstream
pi0 50 2.0 π⁰, 50 GeV, 2 X₀ upstream
pi- 50 2.0 π⁻, 50 GeV, 2 X₀ upstream
e- 50 6.0 e⁻, 50 GeV, 6 X₀ upstream
LIST
rm -rf "$TMP"
