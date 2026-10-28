#!/usr/bin/env python3
"""Build NPxIBD resource website data files from verified sources (no new statistics).
Sources: evidence_matrix.csv (590 genes), pilot shared files (program sets),
SCENIC regulons (scenic_reg.gmt, scenic_regulon_overlap.csv).
Outputs: website/data/*. Verifies counts against manuscript numbers; aborts on mismatch."""
import csv, json, hashlib, shutil, os, sys

BASE = "/home/hatch/workspace/sci-np-ibd"
OUT = os.path.join(BASE, "analysis/09_revision_v1/website/data")
os.makedirs(OUT, exist_ok=True)

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

# 1. evidence matrix: byte-identical copy
src = os.path.join(BASE, "analysis/09_revision_v1/evidence_matrix.csv")
dst = os.path.join(OUT, "evidence_matrix.csv")
shutil.copyfile(src, dst)
with open(dst) as f:
    rows = list(csv.DictReader(f))
assert len(rows) == 590, len(rows)
assert len({r["SYMBOL"] for r in rows}) == 590

# 2. program sets from pilot shared files
def load_shared(p):
    up, dn = [], []
    with open(p) as f:
        for r in csv.DictReader(f):
            (up if int(r["dir_NP"]) == 1 else dn).append(r["SYMBOL"])
    return up, dn

uc_up, uc_dn = load_shared(os.path.join(BASE, "pilot/shared_NPxUC.csv"))
cd_up, cd_dn = load_shared(os.path.join(BASE, "pilot/shared_NPxCD.csv"))
assert (len(uc_up), len(uc_dn)) == (473, 78), (len(uc_up), len(uc_dn))
assert (len(cd_up), len(cd_dn)) == (174, 52), (len(cd_up), len(cd_dn))
core_up = sorted(set(uc_up) & set(cd_up))
core_dn = sorted(set(uc_dn) & set(cd_dn))
assert len(core_up) + len(core_dn) == 187, (len(core_up), len(core_dn))
programs = {"NPxUC_UP": uc_up, "NPxUC_DOWN": uc_dn, "NPxCD_UP": cd_up, "NPxCD_DOWN": cd_dn,
            "core_UP": core_up, "core_DOWN": core_dn}
with open(os.path.join(OUT, "programs.json"), "w") as f:
    json.dump(programs, f)

# 3. SCENIC regulons (GMT -> JSON)
regulons = {}
with open(os.path.join(BASE, "analysis/07_scenic_traj/scenic_reg.gmt")) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        tf = parts[0].replace("_regulon", "")
        regulons[tf] = [g for g in parts[2:] if g]
assert len(regulons) == 278, len(regulons)
with open(os.path.join(OUT, "regulons.json"), "w") as f:
    json.dump(regulons, f)

# 4. regulon-program overlap
overlap = []
with open(os.path.join(BASE, "analysis/07_scenic_traj/scenic_regulon_overlap.csv")) as f:
    for r in csv.DictReader(f):
        overlap.append({"TF": r["TF"], "program": r["program"], "n_target": int(r["n_target"]),
                        "n_overlap": int(r["n_overlap"]), "p": float(r["p"]), "padj": float(r["padj"])})
with open(os.path.join(OUT, "regulon_overlap.json"), "w") as f:
    json.dump(overlap, f)

# 5. downloads dir: build-time copies of the canonical R function + matrix docs
# (canonical source stays analysis/09_revision_v1/; these copies are regenerated here)
DL = os.path.join(BASE, "analysis/09_revision_v1/website/downloads")
os.makedirs(DL, exist_ok=True)
for fn in ["score_program.R", "evidence_matrix_README.md", "evidence_matrix_qc.txt"]:
    shutil.copyfile(os.path.join(BASE, "analysis/09_revision_v1", fn), os.path.join(DL, fn))

with open(os.path.join(OUT, "data_manifest.txt"), "w") as f:
    for fn in ["evidence_matrix.csv", "programs.json", "regulons.json", "regulon_overlap.json"]:
        p = os.path.join(OUT, fn)
        f.write(f"{fn}\t{os.path.getsize(p)}\t{md5(p)}\n")
    f.write(f"counts: matrix=590 UC_UP=473 UC_DOWN=78 CD_UP=174 CD_DOWN=52 core={len(core_up)+len(core_dn)} regulons=278 overlap_rows={len(overlap)}\n")
print("site data OK: matrix=590 UC 473/78 CD 174/52 core=%d regulons=278 overlap_rows=%d" % (len(core_up)+len(core_dn), len(overlap)))
