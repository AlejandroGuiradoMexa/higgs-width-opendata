#!/bin/bash
# Run the g4calo samples used by the shower-shape study, 8 jobs in parallel.
# Usage: ./run_production.sh [n_events]
set -e
N=${1:-1000}
cd "$(dirname "$0")"
source "$HOME/opt/geant4/bin/geant4.sh"
mkdir -p out
{
  for p in e- gamma pi0 pi-; do
    for E in 10 25 50 100; do echo "$p $E 2.0"; done
  done
  for p in e- gamma; do
    for X in 0.0 1.0 4.0 6.0; do echo "$p 50 $X"; done
  done
} | awk -v n="$N" '{seed=1000+NR; printf "%s %s %s %s out/%s_%sGeV_X0%s.txt %d\n", $1, $2, n, $3, $1, $2, $3, seed}' \
  | xargs -P 8 -L 1 sh -c 'test -s "$5" || ./build/calo "$1" "$2" "$3" "$4" "$5" "$6" > /dev/null 2>&1' sh
echo done
