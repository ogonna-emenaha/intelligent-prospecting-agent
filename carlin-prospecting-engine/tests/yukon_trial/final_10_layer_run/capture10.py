import sys, numpy as np, rasterio, geopandas as gpd, json
sys.path.insert(0, "/home/claude/engine/scripts")
from engine import fit_pu_rf, score
from io_utils import read_layers
from layer_registry import IDS
from scipy.cluster.hierarchy import fcluster, linkage
v = gpd.read_file("prep/validation_carlin_ygs.gpkg"); xy = np.c_[v.geometry.x, v.geometry.y]
clu = fcluster(linkage(xy, "single"), 15000, "distance")
L, prof = read_layers("v3/layers", IDS); have = [s for s in IDS if s in L]
tract = np.all(np.isfinite(np.stack([L[s] for s in have])), 0); X = np.stack([L[s][tract] for s in have], 1)
r, c = rasterio.transform.rowcol(prof["transform"], xy[:, 0], xy[:, 1]); r = np.array(r); c = np.array(c)
Xd = np.stack([L[s][r, c] for s in have], 1)
cuts = [0.99, 0.95, 0.90, 0.80]
def table(p, label):
    row = {f"top {round((1-q)*100)}%": f"{int((p >= q).sum())}/{len(p)}" for q in cuts}
    cl = {k: float(np.median(p[clu == k])) for k in np.unique(clu)}
    row.update({f"areas top {round((1-q)*100)}%": f"{sum(x >= q for x in cl.values())}/{len(cl)}" for q in cuts})
    print(label, json.dumps(row)); return row
res = {}
# 1. zero-shot Nevada archetype (run R5, same 8 layers)
with rasterio.open("v3/zero_shot_magmatic_10/carlin_similarity_percentile.tif") as ds: z = ds.read(1)
res["zero_shot_nevada"] = table(z[r, c], "zero-shot Nevada")
# 2. local Yukon model, leave one area out (honest check)
rng = np.random.default_rng(0); bgi = rng.choice(len(X), 20000, replace=False); sub = rng.choice(len(X), 20000, replace=False)
p = np.full(len(xy), np.nan)
for k in np.unique(clu):
    te = clu == k; m = fit_pu_rf(Xd[~te], X[bgi], n_bags=25, seed=int(k))
    sr = score(m, X[sub]); st = score(m, Xd[te]); p[te] = [(sr < s).mean() for s in st]
res["local_yukon_leave_one_area_out"] = table(p, "local Yukon, held-out areas")
# 3. local Yukon model, full fit (deposits scored by a model that saw them; for reference only)
with rasterio.open("v3/local_yukon_10/local_rf_percentile.tif") as ds: lf = ds.read(1)
res["local_yukon_full_fit_not_a_test"] = table(lf[r, c], "local Yukon, trained on all")
import pandas as pd
d = pd.DataFrame({"name": v.NAME, "area": clu, "zero_shot": z[r, c].round(3), "local_held_out": p.round(3), "local_full_fit": lf[r, c].round(3)})
print(d.groupby("area").agg(n=("name","size"), names=("name", lambda s: ", ".join(s[:3])), zero_shot=("zero_shot","median"), local_held_out=("local_held_out","median")).round(2).to_string())
d.to_csv("v3/capture_by_deposit_10.csv", index=False); json.dump(res, open("v3/capture_summary_10.json", "w"), indent=1)
