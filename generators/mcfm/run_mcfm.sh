#!/bin/bash
# Run all MCFM cards in parallel inside the MCFM Bin directory (it needs Pdfdata/ there).
# Cards are copied with short names because MCFM truncates long input paths.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HOME/opt/src/MCFM-7.0_JHUGen/Bin"
mkdir -p "$HERE/output"
for card in "$HERE"/*.DAT; do
  name=$(basename "$card" .DAT)
  cp "$card" "$BIN/$name.DAT"
  ( cd "$BIN" && ./mcfm "$name.DAT" > "$HERE/output/$name.log" 2>&1 && mv ./*_"$name".dat "$HERE/output/" ) &
done
wait
echo done
