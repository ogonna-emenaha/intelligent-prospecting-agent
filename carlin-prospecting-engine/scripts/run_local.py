#!/usr/bin/env python3
"""Mode B: train a local prospectivity model from the user's own deposits on the user's own evidence layers.
Usage: python run_local.py --layers DIR --deposits points.gpkg --out DIR [--model rf|ocsvm|xgboost] [--compare-nevada]"""
import os, sys, json, argparse, warnings; warnings.filterwarnings("ignore")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from engine import fit_pu_rf, score, percentile_rank
from layer_registry import IDS, NAME
from io_utils import read_layers, write_raster, deposits_xy

def run(layers_dir, deposits_xy, out_dir, model="rf", n_bags=25, min_deposits=10, seed=0, compare_nevada=False, archetype="magmatic", tract_basis=None, allow_unfiltered=False):
    import rasterio
    from sklearn.cluster import KMeans
    os.makedirs(out_dir, exist_ok=True)
    L, prof = read_layers(layers_dir, IDS); have = [s for s in IDS if s in L]
    tract = np.all(np.isfinite(np.stack([L[s] for s in have])), 0); X = np.stack([L[s][tract] for s in have], 1)
    r, c = rasterio.transform.rowcol(prof["transform"], deposits_xy[:, 0], deposits_xy[:, 1]); r = np.array(r); c = np.array(c)
    ok = (r >= 0) & (r < tract.shape[0]) & (c >= 0) & (c < tract.shape[1]); ok[ok] = tract[r[ok], c[ok]]
    Xd = np.stack([L[s][r[ok], c[ok]] for s in have], 1); n = len(Xd)
    rep = dict(layers_used=[NAME[s] for s in have], n_deposits_in_area=int(n))
    if n < min_deposits:
        rep["verdict"] = f"Too few deposits ({n}) for a reliable local model; minimum {min_deposits} (provisional). Use zero-shot mode, or add deposits."
        json.dump(rep, open(os.path.join(out_dir, "local_summary.json"), "w"), indent=1); return rep
    rng = np.random.default_rng(seed); bgi = rng.choice(len(X), min(20000, len(X)), replace=False)
    def fit_predict(Xtr, Xte_list):
        if model == "rf":
            m = fit_pu_rf(Xtr, X[bgi], n_bags=n_bags, seed=seed); return [score(m, t) for t in Xte_list]
        if model == "ocsvm":
            from sklearn.svm import OneClassSVM
            m = OneClassSVM(kernel="rbf", nu=0.15, gamma="scale").fit(Xtr); return [m.decision_function(t) for t in Xte_list]
        if model == "xgboost":
            from xgboost import XGBClassifier
            ni = rng.choice(bgi, len(Xtr), replace=False); m = XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.1, verbosity=0).fit(np.vstack([Xtr, X[ni]]), np.r_[np.ones(len(Xtr)), np.zeros(len(Xtr))]); return [m.predict_proba(t)[:, 1] for t in Xte_list]
        raise ValueError(model)
    prob, = fit_predict(Xd, [X]); pr = percentile_rank(prob)
    g = np.full(tract.shape, np.nan); g[tract] = prob; write_raster(os.path.join(out_dir, f"local_{model}_probability.tif"), g, prof)
    g2 = np.full(tract.shape, np.nan); g2[tract] = pr; write_raster(os.path.join(out_dir, f"local_{model}_percentile.tif"), g2, prof)
    # spatial-block cross-validation of the user's own deposits (held-out blocks)
    if n >= 10:
        k = min(5, n // 3); folds = KMeans(k, random_state=0, n_init=10).fit_predict(deposits_xy[ok]); pcts = []
        sub = rng.choice(len(X), min(20000, len(X)), replace=False)
        for f in range(k):
            tr = folds != f; te = folds == f
            s_ref, s_te = fit_predict(Xd[tr], [X[sub], Xd[te]]); pcts += [float((s_ref < v).mean()) for v in s_te]
        p = np.array(pcts); rep["spatial_block_cv"] = dict(folds=k, mean_percentile_of_held_out_deposits=round(float(p.mean()), 3), share_in_top_20pct=round(float((p >= 0.8).mean()), 3), chance_share=0.2)
    if compare_nevada:
        import run_zero_shot
        z = run_zero_shot.run(layers_dir, os.path.join(out_dir, "zero_shot"), deposits_xy=deposits_xy, archetype=archetype, tract_basis=tract_basis, allow_unfiltered=allow_unfiltered)
        if z.get("verdict", "").startswith("Not run"): rep["agreement_with_nevada_signature"] = dict(note=z["reason"]); rep["verdict"] = "Local model trained."; json.dump(rep, open(os.path.join(out_dir, "local_summary.json"), "w"), indent=1); return rep
        import rasterio as rio
        with rio.open(os.path.join(out_dir, "zero_shot", "carlin_similarity_percentile.tif")) as ds: zs = ds.read(1)
        zv = zs[tract]; from scipy.stats import spearmanr
        rep["agreement_with_nevada_signature"] = dict(spearman=round(float(spearmanr(pr, zv).correlation), 3),
            both_top20=round(float(((pr >= .8) & (zv >= .8)).sum() / max((pr >= .8).sum(), 1)), 3))
    rep["verdict"] = "Local model trained."; json.dump(rep, open(os.path.join(out_dir, "local_summary.json"), "w"), indent=1); return rep

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--layers", required=True); ap.add_argument("--deposits", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--model", default="rf", choices=["rf", "ocsvm", "xgboost"]); ap.add_argument("--min-deposits", type=int, default=10); ap.add_argument("--compare-nevada", action="store_true"); ap.add_argument("--archetype", default="magmatic", choices=["magmatic", "magma_poor", "older_buried"]); ap.add_argument("--tract-basis"); ap.add_argument("--allow-unfiltered", action="store_true")
    a = ap.parse_args()
    xy = deposits_xy(a.deposits, a.layers)
    print(json.dumps(run(a.layers, xy, a.out, a.model, min_deposits=a.min_deposits, compare_nevada=a.compare_nevada, archetype=a.archetype, tract_basis=a.tract_basis, allow_unfiltered=a.allow_unfiltered), indent=1))
