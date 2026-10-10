#!/usr/bin/env python3
"""Stricter evaluation protocol for a Carlin run in a district that has known deposits. Use it to find out whether the Nevada archetype adds anything,
and whether a good-looking score is better than luck. It never changes the map; it only tests.

What it does (all on the user's layers, with the user's deposits grouped into AREAS by single-linkage distance):
  1. Baseline        a deliberately simple map (reactive host rock + fault density, equal rank weights). Any Carlin model has to beat this.
  2. Zero-shot       each chosen Nevada archetype applied with no local deposits (k = 0).
  3. Local only      a model trained on k other areas' deposits only (k = 1, 2, ... and all other areas), scored on the held-out area.
  4. Nevada + local  the archetype's deposits PLUS the same k other areas, scored on the held-out area. If this beats 'local only' at the same k, the Nevada
                     archetype transfers value when local data is scarce, even if zero-shot is weak.
  5. Negative control  deposits of other types (for example MVT, SEDEX, orogenic gold) scored by the same maps; AUC of Carlin against each other-type group.
Whole areas are held out, never single points. Percentiles are measured against a random sample of the study area's cells.
Usage: python run_evaluation.py --layers DIR --positives carlin.gpkg --out DIR --tag zeroshot --archetypes magmatic,magma_poor [--negatives other.gpkg --negative-group-field GROUP]
       python run_evaluation.py --layers DIR --positives carlin.gpkg --out DIR --tag fs_magma_poor --archetypes magma_poor --fewshot 1,2
       python run_evaluation.py --report --out DIR            (merges every eval_*.json in DIR)"""
import os, sys, json, glob, argparse, itertools, warnings; warnings.filterwarnings("ignore")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from engine import fit_pu_rf, score
from layer_registry import IDS, NAME, DIAGNOSTIC
from io_utils import read_layers, deposits_xy
from run_zero_shot import load_archetype, REF
CUTS = (0.99, 0.95, 0.90, 0.80)

def pct(vals, ref_sorted): return np.searchsorted(ref_sorted, vals, side="left") / len(ref_sorted)
def metrics(p_by_area):
    """p_by_area: dict area -> array of deposit percentiles. Returns deposit-level and area-level summaries."""
    allp = np.concatenate(list(p_by_area.values())); med = {a: float(np.median(v)) for a, v in p_by_area.items()}
    out = {f"deposits_top{round((1-c)*100)}pct": f"{int((allp >= c).sum())}/{len(allp)}" for c in CUTS}
    out.update({f"areas_top{round((1-c)*100)}pct": f"{sum(m >= c for m in med.values())}/{len(med)}" for c in CUTS})
    out.update(mean_percentile_deposits=round(float(allp.mean()), 3), mean_percentile_areas=round(float(np.mean(list(med.values()))), 3), n_deposits=int(len(allp)), n_areas=len(med))
    return out
def auc(a, b):
    a = np.asarray(a); b = np.asarray(b); return float(((a[:, None] > b[None, :]).sum() + 0.5 * (a[:, None] == b[None, :]).sum()) / (len(a) * len(b)))

def run(layers_dir, positives, out_dir, tag, archetypes=("magmatic",), fewshot=(), negatives=None, neg_field=None, area_km=15.0, exclude=(), n_cells=30000, seed=0, max_combos=6, local_share=0.5, bags=8, trees=40, baseline=True):
    import rasterio, geopandas as gpd
    from scipy.cluster.hierarchy import fcluster, linkage
    os.makedirs(out_dir, exist_ok=True); rng = np.random.default_rng(seed)
    L, prof = read_layers(layers_dir, IDS); allhave = [s for s in IDS if s in L and s not in exclude]
    tract = np.all(np.isfinite(np.stack([L[s] for s in allhave])), 0)
    def cells_of(xy):
        r, c = rasterio.transform.rowcol(prof["transform"], xy[:, 0], xy[:, 1]); r = np.array(r); c = np.array(c); ok = (r >= 0) & (r < tract.shape[0]) & (c >= 0) & (c < tract.shape[1]); ok[ok] = tract[r[ok], c[ok]]; return r, c, ok
    xy = deposits_xy(positives, layers_dir); r, c, ok = cells_of(xy); xy = xy[ok]; r, c = r[ok], c[ok]
    area = fcluster(linkage(xy, "single"), area_km * 1000, "distance"); areas = np.unique(area)
    idx = np.flatnonzero(tract.ravel()); samp = rng.choice(idx, min(n_cells, len(idx)), replace=False); sr, sc = np.unravel_index(samp, tract.shape)
    def feat(cols, rr, cc): return np.stack([L[s][rr, cc] for s in cols], 1)
    res = dict(tag=tag, n_deposits=int(len(xy)), n_deposits_dropped_outside_tract=int((~ok).sum()), areas={int(a): int((area == a).sum()) for a in areas}, area_km=area_km, cells_in_study_area=int(tract.sum()), conditions={}, negative_control={})
    neg = None
    if negatives:
        ng = gpd.read_file(negatives); crs = rasterio.open(sorted(glob.glob(os.path.join(layers_dir, "*.tif")))[0]).crs; ng = ng.to_crs(crs)
        nxy = np.c_[ng.geometry.x.values, ng.geometry.y.values]; nr, nc, nok = cells_of(nxy); grp = ng[neg_field].astype(str).values[nok] if neg_field else np.array(["other"] * int(nok.sum())); neg = (nr[nok], nc[nok], grp)
        res["negatives_in_study_area"] = {g: int((grp == g).sum()) for g in np.unique(grp)}
    def neg_report(name, scorer):
        if neg is None: return
        pr_ref = np.sort(scorer(sr, sc)); pdep = pct(scorer(r, c), pr_ref); pn = pct(scorer(neg[0], neg[1]), pr_ref)
        res["negative_control"][name] = {g: dict(n=int((neg[2] == g).sum()), carlin_vs_group_AUC=round(auc(pdep, pn[neg[2] == g]), 3), group_mean_percentile=round(float(pn[neg[2] == g].mean()), 3)) for g in np.unique(neg[2])}
        res["negative_control"][name]["_carlin_mean_percentile"] = round(float(pdep.mean()), 3)
    # 1. baseline
    if baseline and "trap_reactive_carbonate" in L and "arch_fault_density" in L:
        def base(rr, cc): return 0.5 * L["trap_reactive_carbonate"][rr, cc] + 0.5 * L["arch_fault_density"][rr, cc]
        ref = np.sort(base(sr, sc)); pdep = pct(base(r, c), ref); res["conditions"]["baseline: reactive host + fault density"] = metrics({int(a): pdep[area == a] for a in areas}); neg_report("baseline: reactive host + fault density", base)
    # 2-4 archetypes
    for a_name in archetypes:
        am, Xn, _, aids, _ = load_archetype(a_name); skip = set(exclude) | set(am["exclude_layers"]); have = [s for s in IDS if s in L and s not in skip and s in aids]
        ncol = [aids.index(s) for s in have]; Xnev = Xn[:, ncol]; bgz = np.load(os.path.join(REF, "reference_background.npz"), allow_pickle=True); Xbg_n = bgz["X"][:, [list(map(str, bgz["ids"])).index(s) for s in have]].astype(float)
        Xs = feat(have, sr, sc); Xd = feat(have, r, c); Xn_pts = feat(have, neg[0], neg[1]) if neg else None
        res.setdefault("layers_used", {})[a_name] = [NAME[s] for s in have]; res.setdefault("diagnostic_layers_used", {})[a_name] = [NAME[s] for s in have if s in DIAGNOSTIC]
        # layer sign check: does each layer separate deposits from background the same way in this district as in the Nevada archetype?
        from sklearn.metrics import roc_auc_score
        sc_rows = []
        for j, sname in enumerate(have):
            ad = roc_auc_score(np.r_[np.ones(len(Xd)), np.zeros(len(Xs))], np.r_[Xd[:, j], Xs[:, j]]); an = roc_auc_score(np.r_[np.ones(len(Xnev)), np.zeros(len(Xbg_n))], np.r_[Xnev[:, j], Xbg_n[:, j]])
            flag = "SIGN CONFLICT: check how this layer was built, or leave it out" if (ad - 0.5) * (an - 0.5) < 0 and max(abs(ad - 0.5), abs(an - 0.5)) > 0.15 and min(abs(ad - 0.5), abs(an - 0.5)) > 0.05 else ("weak here" if abs(ad - 0.5) < 0.1 and abs(an - 0.5) > 0.2 else "")
            sc_rows.append(dict(layer=NAME[sname], district_AUC=round(float(ad), 2), archetype_AUC=round(float(an), 2), flag=flag))
        res.setdefault("layer_sign_check", {})[a_name] = sc_rows
        # zero-shot
        m0 = fit_pu_rf(Xnev, Xbg_n, n_bags=25, seed=seed); ref = np.sort(score(m0, Xs)); pdep = pct(score(m0, Xd), ref)
        res["conditions"][f"zero-shot {a_name} (k=0)"] = metrics({int(a): pdep[area == a] for a in areas}); neg_report(f"zero-shot {a_name}", lambda rr, cc, m=m0, cols=have: score(m, feat(cols, rr, cc)))
        # few-shot: local only vs Nevada + local, whole areas held out
        for k in fewshot:
            kk = len(areas) - 1 if k == "all" else int(k)
            if kk > len(areas) - 1: continue
            for kind in (f"local only [{len(have)} layers of {a_name} set]", f"Nevada({a_name}) + local"):
                pa = {}
                for h in areas:
                    others = [x for x in areas if x != h]; combos = list(itertools.combinations(others, kk))
                    if len(combos) > max_combos: combos = [combos[i] for i in rng.choice(len(combos), max_combos, replace=False)]
                    acc = np.zeros((int((area == h).sum()), 0)); ps = []
                    for cb in combos:
                        loc = np.isin(area, cb); Xl = Xd[loc]; Xcell = Xs
                        if kind.startswith("local only"): m = fit_pu_rf(Xl, Xcell, n_bags=bags, n_trees=trees, seed=seed)
                        else:
                            rep = int(np.clip(np.ceil(len(Xnev) * local_share / (1 - local_share) / len(Xl)), 1, 20)); Xp = np.vstack([Xnev] + [Xl] * rep)
                            m = fit_pu_rf(Xp, np.vstack([Xbg_n[rng.choice(len(Xbg_n), min(10000, len(Xbg_n)), replace=False)], Xcell]), n_bags=bags, n_trees=trees, seed=seed)
                        ref = np.sort(score(m, Xcell)); ps.append(pct(score(m, Xd[area == h]), ref))
                    pa[int(h)] = np.mean(ps, axis=0)
                res["conditions"][f"{kind} (k={kk} other areas)"] = metrics(pa)
    json.dump(res, open(os.path.join(out_dir, f"eval_{tag}.json"), "w"), indent=1); return res

def report(out_dir):
    import pandas as pd
    rows = []; neg = []
    for f in sorted(glob.glob(os.path.join(out_dir, "eval_*.json"))):
        j = json.load(open(f))
        for k, v in j["conditions"].items(): rows.append(dict(source=j["tag"], condition=k, **v))
        for k, v in j["negative_control"].items():
            for g, s in v.items():
                if not g.startswith("_"): neg.append(dict(source=j["tag"], model=k, group=g, **s))
    d = pd.DataFrame(rows).drop_duplicates("condition"); d.to_csv(os.path.join(out_dir, "evaluation_conditions.csv"), index=False); pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); print(d.drop(columns=["source", "n_areas", "n_deposits"]).to_string(index=False))
    if neg: n = pd.DataFrame(neg).drop_duplicates(["model", "group"]); n.to_csv(os.path.join(out_dir, "evaluation_negative_control.csv"), index=False); print(); print(n.drop(columns="source").to_string(index=False))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--layers"); ap.add_argument("--positives"); ap.add_argument("--out", required=True); ap.add_argument("--tag", default="run")
    ap.add_argument("--archetypes", default="magmatic"); ap.add_argument("--fewshot", default="", help="comma list of k values, 'all' allowed, e.g. 1,2,all"); ap.add_argument("--negatives"); ap.add_argument("--negative-group-field")
    ap.add_argument("--area-km", type=float, default=15.0); ap.add_argument("--exclude", default=""); ap.add_argument("--max-combos", type=int, default=6); ap.add_argument("--no-baseline", action="store_true"); ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.report: report(a.out)
    else:
        fs = [x if x == "all" else int(x) for x in a.fewshot.split(",") if x]
        r = run(a.layers, a.positives, a.out, a.tag, a.archetypes.split(","), fs, a.negatives, a.negative_group_field, a.area_km, [x for x in a.exclude.split(",") if x], max_combos=a.max_combos, baseline=not a.no_baseline)
        print(json.dumps({k: v for k, v in r.items() if k in ("n_deposits", "areas", "conditions", "negative_control")}, indent=1))
