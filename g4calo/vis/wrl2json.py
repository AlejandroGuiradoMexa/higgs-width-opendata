"""Convert a Geant4 VRML2FILE output (geometry + trajectories) into the compact JSON read by the 3D view of the
project page (site/geant4.js).

The VRML is written by the VRML 2.0 driver of Geant4 (macros/vrml.mac): one IndexedFaceSet per placed volume and
one IndexedLineSet per trajectory, coloured by particle. Here the volumes are reduced to their bounding boxes and
the trajectories are thinned for the web: each polyline is simplified (Ramer-Douglas-Peucker, 0.5 mm tolerance)
and rounded to 0.5 mm, tracks shorter than `min_len` mm are dropped, and everything before `z_min` is cut (the gun is
1.5 m upstream).

Usage: python wrl2json.py g4_00.wrl out.json [--label "e-, 50 GeV"]
"""
import argparse
import json
import math
import re

COLOURS = {"1 0 0": "e-", "0 0 1": "e+", "0 1 0": "gamma"}  # macros/vrml.mac; anything else is "other"
Q = 0.5  # mm per integer unit in the output


def simplify(pts, tol):
    """Ramer-Douglas-Peucker: drop points closer than `tol` to the chord of their neighbours."""
    if len(pts) < 3:
        return pts
    a, b = pts[0], pts[-1]
    ab = [b[i] - a[i] for i in range(3)]
    nab = math.sqrt(sum(c * c for c in ab)) or 1e-12
    best, idx = -1.0, 0
    for i in range(1, len(pts) - 1):
        ap = [pts[i][k] - a[k] for k in range(3)]
        cr = (ap[1] * ab[2] - ap[2] * ab[1], ap[2] * ab[0] - ap[0] * ab[2], ap[0] * ab[1] - ap[1] * ab[0])
        d = math.sqrt(sum(c * c for c in cr)) / nab
        if d > best:
            best, idx = d, i
    if best <= tol:
        return [a, b]
    return simplify(pts[:idx + 1], tol)[:-1] + simplify(pts[idx:], tol)


def blocks(text):
    """Yield (kind, name, colour, points) for each SOLID / POLYLINE block."""
    for m in re.finditer(r"#---------- (SOLID: (\S+)|POLYLINE)\n(.*?)(?=#----------|#End of file)", text, re.S):
        body = m.group(3)
        col = re.search(r"diffuseColor ([\d.e+-]+ [\d.e+-]+ [\d.e+-]+)", body)
        pts = re.search(r"point \[(.*?)\]", body, re.S)
        xyz = [tuple(map(float, p.split())) for p in pts.group(1).split(",") if p.strip()] if pts else []
        yield ("solid" if m.group(2) else "track"), m.group(2), col.group(1) if col else "", xyz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("wrl")
    ap.add_argument("out")
    ap.add_argument("--label", default="")
    ap.add_argument("--min-len", type=float, default=3.0)
    ap.add_argument("--z-min", type=float, default=-400.0)
    a = ap.parse_args()

    text = open(a.wrl).read()
    solids, tracks, n_in = [], {}, 0
    for kind, name, col, xyz in blocks(text):
        if kind == "solid":
            xs, ys, zs = zip(*xyz)
            solids.append([name.split(":")[0], round(min(xs)), round(max(xs)), round(min(ys)), round(max(ys)),
                           round(min(zs), 1), round(max(zs), 1)])
            continue
        n_in += 1
        pts = [p for p in xyz if p[2] >= a.z_min]
        if len(pts) < 2:
            continue
        length = sum(math.dist(pts[i], pts[i + 1]) for i in range(len(pts) - 1))
        if length < a.min_len:
            continue
        kept = simplify(pts, 0.5)
        flat = [round(c / Q) for p in kept for c in p]
        tracks.setdefault(COLOURS.get(col, "other"), []).append(flat)

    out = {"label": a.label, "units_mm": Q, "solids": solids, "tracks": tracks,
           "n_tracks_in": n_in, "n_tracks_kept": sum(len(v) for v in tracks.values())}
    with open(a.out, "w") as f:
        json.dump(out, f, separators=(",", ":"))
    print(f"{a.out}: {len(solids)} volumes, {out['n_tracks_kept']} of {n_in} tracks kept")


if __name__ == "__main__":
    main()
