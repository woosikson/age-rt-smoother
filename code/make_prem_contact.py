#!/usr/bin/env python3
# age-rt-smoother / code/make_prem_contact.py — alternative contact matrices (paper Fig S8, S9)
# ─────────────────────────────────────────────────────────────────────────────
# Prem 2017 / 2021 synthetic contact matrices for Korea (16x16, five-year bands)
#   -> the six bands used here (0-5 / 6-11 / 12-17 / 18-44 / 45-64 / 65+).
#   (1) aggregation: weighted by the single-year population of Pops_Dec2022.csv so that the
#       six-band boundaries are exact. The only approximation is that contacts are assumed
#       uniform within a five-year band.
#   (2) transposition: Prem reports row = ego (participant); the convention here is
#       column = infector, so C_prem = M_agg^T.
# Output: data/contact_prem{2017,2021}.csv (6x6, no header)
# The original releases are not redistributed here — download them from the two papers and
# pass their directories:
#   python3 code/make_prem_contact.py --prem2017 <dir with MUestimates_*.xlsx> \
#                                     --prem2021 <dir with contact_all.rdata>
#   (defaults: external/Prem_2017, external/Prem_2021; PREM2017_DIR / PREM2021_DIR also work)
# Prem 2021 ships an R .rdata file, so Rscript is called once to export the KOR matrix as CSV.
# ─────────────────────────────────────────────────────────────────────────────
import csv, io, subprocess, os, numpy as np, pandas as pd
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
DDIR = os.path.join(ROOT, "data")
PREM = {"--prem2017": os.environ.get("PREM2017_DIR", os.path.join(ROOT, "external", "Prem_2017")),
        "--prem2021": os.environ.get("PREM2021_DIR", os.path.join(ROOT, "external", "Prem_2021"))}
_a = sys.argv[1:]
while _a:
    if _a[0] not in PREM or len(_a) < 2:
        raise SystemExit(f"usage: make_prem_contact.py [--prem2017 DIR] [--prem2021 DIR]")
    PREM[_a[0]] = _a[1]; _a = _a[2:]
for k, v in PREM.items():
    if not os.path.isdir(v):
        raise SystemExit(f"{k} directory not found: {v}\n  (see data/README.md for where to download it)")

# --- single-year population, national row of Pops_Dec2022.csv (EUC-KR encoded) ---
# The published file has Korean column headers. They are written as \u escapes so that this
# source stays pure ASCII:
#   \uc804\uad6d          "national total" (the row we take)
#   2022\ub144_\uacc4_    "2022_total_" prefix of every single-year population column
#   \uc5f0\ub839\uad6c\uac04 / \ucd1d\uc778\uad6c   "age range" / "total population" columns, both skipped
#   \uc138 / \uc774\uc0c1   "years old" / "and over" suffixes
rows = list(csv.reader(io.open(os.path.join(DDIR, "Pops_Dec2022.csv"), encoding="euc-kr")))
hdr = rows[0]; nat = [r for r in rows if r[0].startswith("\uc804\uad6d")][0]
num = lambda x: int(x.replace(",", ""))
Nyear = {}
for i, h in enumerate(hdr):
    if h.startswith("\u0032\u0030\u0032\u0032\ub144\u005f\uacc4\u005f") and "\uc5f0\ub839\uad6c\uac04" not in h and "\ucd1d\uc778\uad6c" not in h:
        lab = h.replace("\u0032\u0030\u0032\u0032\ub144\u005f\uacc4\u005f", "")
        if lab.endswith("\uc138") and "\uc774\uc0c1" not in lab: a = int(lab.replace("\uc138", ""))
        elif "\u0031\u0030\u0030\uc138\u0020\uc774\uc0c1" in lab: a = 100
        else: continue
        Nyear[a] = num(nat[i])
AGES = sorted(Nyear)
band = lambda a: min(a // 5, 15)                       # Prem's 16 bands (0-4 .. 75-79); 80+ folded into 75-79
grp = lambda a: 0 if a <= 5 else 1 if a <= 11 else 2 if a <= 17 else 3 if a <= 44 else 4 if a <= 64 else 5
Nband = {}
for a in AGES: Nband[band(a)] = Nband.get(band(a), 0) + Nyear[a]

def aggregate_transpose(M):                            # M: 16x16 Prem (row = ego) -> 6x6, column = infector
    Magg = np.zeros((6, 6))
    for I in range(6):
        egos = [y for y in AGES if grp(y) == I]; den = sum(Nyear[y] for y in egos)
        for J in range(6):
            tgt = [z for z in AGES if grp(z) == J]
            s = sum(Nyear[y] * M[band(y), band(z)] * Nyear[z] / Nband[band(z)] for y in egos for z in tgt)
            Magg[I, J] = s / den
    return Magg.T                                      # transpose → col=ego(infector)

# --- Prem 2017 (xlsx) ---
M17 = pd.read_excel(os.path.join(PREM["--prem2017"], "MUestimates_all_locations_2.xlsx"),
                    sheet_name="Republic of Korea", header=None).values.astype(float)
# --- Prem 2021 (rdata → csv) ---
kor21 = os.path.join(DDIR, "_prem2021_kor.csv")            # cached intermediate, reused on re-run
if not os.path.exists(kor21):
    rdata = os.path.join(PREM["--prem2021"], "output", "syntheticmatrices", "contact_all.rdata")
    if not os.path.exists(rdata):
        rdata = os.path.join(PREM["--prem2021"], "contact_all.rdata")
    subprocess.run(["Rscript", "-e",
        f'load("{rdata}"); write.csv(contact_all[["KOR"]], "{kor21}", row.names=FALSE)'], check=True)
M21 = pd.read_csv(kor21).values.astype(float)

for name, M in [("contact_prem2017", M17), ("contact_prem2021", M21)]:
    C = aggregate_transpose(M)
    pd.DataFrame(C).to_csv(os.path.join(DDIR, f"{name}.csv"), index=False, header=False)
    print(f"saved data/{name}.csv  (diagonal self-contacts: {np.round(np.diag(C), 2)})")
