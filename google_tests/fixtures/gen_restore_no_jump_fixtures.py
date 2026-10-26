#!/usr/bin/env python3
"""Regenerate the restoreNoJump fixtures at full precision.

Source: res/10nsstart+fitting/10ns.trr -- the 10 ns frame of the ubiquitin run, which is where the
committed 2-decimal-place fixtures came from (established by Kabsch-fitting the committed
ensemble_coord_10ns.csv onto 10ns.pdb atoms 1..1231: RMSD 0.0050 A, i.e. the rounding floor, against
2.35 A / 2.29 A for the 00ns and 02ns frames). The .trr is used rather than the .gro or .pdb because
it carries full float32 precision, where the .gro rounds to 0.001 nm and the .pdb to 0.001 A.

Output:
  ensemble_coord_10ns.csv                 atoms 1..1250 of that frame (protein 1..1231 + 19 waters),
                                          in Angstrom, shortest round-tripping decimal.
  ensemble_coord_moldel2_PBC_Broken.csv   the same, with the SAME 80 z-coordinates as the committed
                                          fixture moved down by exactly one box vector in z.

The frame is used as it comes out of the trajectory, NOT in the rigidly-fitted orientation of the
committed fixture: the fit cannot be recovered better than the 2-dp data it would have to be
estimated from, whereas the raw frame is exactly reproducible from a tracked file. It is also the
frame in which a periodic wrap is physically meaningful, the fitted one having been rotated out of
the box frame. Atom indexing and order are unchanged, so the test's guide/sub index lists still
name the same atoms.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trr import read_trr

repo = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
NATOMS = 1250

r = read_trr(os.path.join(repo, "res/10nsstart+fitting/10ns.trr"))
assert r["natoms"] == 25951, r["natoms"]
box = r["box"]
assert box[0][1] == box[0][2] == box[1][0] == box[1][2] == box[2][0] == box[2][1] == 0.0, "expected a rectangular box"
box_nm_z = float(box[2][2])
box_A_z  = box_nm_z * 10.0
X = r["x"][:NATOMS] * 10.0          # nm -> Angstrom

# The 80 z-coordinates to tear, 0-based rows. Recovered once from the original two-decimal-place
# fixtures (the rows where they differed) and written out here so this script does not depend on
# the fixtures it generates and stays reproducible from tracked files alone. All 80 are inside the
# protein, atoms 1..1231.
TORN = [
    500, 585, 586, 587, 588, 589, 590, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 628,
    629, 630, 631, 632, 633, 634, 635, 1168, 1169, 1170, 1177, 1178, 1179, 1180, 1181, 1182,
    1183, 1184, 1185, 1186, 1187, 1188, 1189, 1192, 1193, 1194, 1195, 1196, 1197, 1198, 1199,
    1200, 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1210, 1211, 1212, 1213, 1214,
    1215, 1216, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225, 1226, 1227, 1228, 1229,
    1230,
]
torn = TORN
assert len(torn) == 80, len(torn)
assert max(torn) < 1231, "the tear must stay inside the protein"

def write(path, rows):
    with open(path, "w") as fh:
        for row in rows:
            fh.write(",".join(repr(float(v)) for v in row) + "\n")

whole = [list(row) for row in X]
write(os.path.join(repo, "res/google_tests/ensemble_coord_10ns.csv"), whole)
broken = [list(row) for row in whole]
for i in torn:
    broken[i][2] = whole[i][2] - box_A_z
write(os.path.join(repo, "res/google_tests/ensemble_coord_moldel2_PBC_Broken.csv"), broken)

print(f"box z            : {box_nm_z!r} nm = {box_A_z!r} A")
print(f"atoms written    : {NATOMS} (protein 1..1231 + waters 1232..1250)")
print(f"z-coords torn    : {len(torn)}  (indices {torn[0]}..{torn[-1]})")
print(f"max |coord|      : {max(abs(v) for row in whole for v in row):.6f} A")
