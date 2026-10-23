#!/bin/bash
# Run the EFT MCFM cards listed in eft/jobs.txt ("[ENV=value ...] name"), 5 at a time,
# inside the MCFM Bin directory.
HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HOME/opt/src/MCFM-7.0_JHUGen/Bin"
OUT="$HERE/eft/output"
mkdir -p "$OUT"
cp "$HERE"/eft/*.DAT "$BIN/"
cd "$BIN"
while read -r line; do
  name="${line##* }"
  envs="${line% *}"
  [ "$envs" = "$line" ] && envs=""
  ( env $envs ./mcfm "$name.DAT" > "$OUT/$name.log" 2>&1 && mv ./*_"$name".dat "$OUT/" ) &
  while [ "$(jobs -rp | wc -l)" -ge 5 ]; do sleep 2; done
done < "$HERE/eft/jobs.txt"
wait
echo done
