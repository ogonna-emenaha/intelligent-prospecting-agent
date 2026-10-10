"""Check that fixed (Nevada-calibrated) scaling reproduces the area-relative Nevada layers exactly on Nevada, so the stored reference stays valid.
Needs the Terra-X Nevada project files (FA). Run from the skill folder: python tests/test_fixed_scale_nevada.py"""
import os, sys, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, rasterio
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "scripts")); sys.path.insert(0, os.path.join(HERE, "..", "tools"))
import build_layers, calibrate_nevada_scale as C
OUT = os.environ.get("FIX_OUT", "/tmp/fix_nevada")
sc = json.load(open(os.path.join(HERE, "..", "reference", "nevada_scale.json"))); ids = list(sc)
res = {}
for mode in ("relative", "fixed"):
    cfg = C.nevada_cfg(mode); cfg["layers"] = {k: v for k, v in cfg["layers"].items() if k in ids}; log = build_layers.build(cfg, f"{OUT}/{mode}"); res[mode] = log
worst = 0
for lid in ids:
    a = rasterio.open(f"{OUT}/relative/{lid}.tif").read(1); b = rasterio.open(f"{OUT}/fixed/{lid}.tif").read(1); m = (a > -9000) & (b > -9000) & np.isfinite(a) & np.isfinite(b)
    d = float(np.abs(a[m] - b[m]).max()); worst = max(worst, d); print(f"{lid:46s} max |fixed - relative| = {d:.2e}   scaling: {res['fixed'][lid]['scaling'][:5]}")
print("WORST", worst); assert worst < 1e-4, "fixed scaling does not reproduce the Nevada layers"
