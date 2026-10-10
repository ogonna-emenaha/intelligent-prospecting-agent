#!/usr/bin/env python3
"""Build the Nevada reference archetypes from the 84 Nevada deposits and test them inside Nevada.
Archetypes (clusters from the original Terra-X clustering): magmatic (54-member cluster, 40 pass the falsification envelope), magma_poor (18), older_buried (12).
Writes reference/reference_archetypes.npz, reference/reference_model_<archetype>.joblib, reference/archetypes.json and tests/archetype_cv_results.csv."""
import os, sys, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, joblib, geopandas as gpd
from sklearn.cluster import KMeans
from sklearn.metrics import roc_auc_score
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, ".."); sys.path.insert(0, os.path.join(ROOT, "scripts"))
from engine import fit_pu_rf, score
from layer_registry import IDS, NAME
FA = os.environ.get("FA", "."); OLD = f"{FA}/Analytical Results/Old Results"
REF = os.path.join(ROOT, "reference")
d = np.load(f"{REF}/reference_deposits.npz", allow_pickle=True); bg = np.load(f"{REF}/reference_background.npz", allow_pickle=True)
X84 = d["X_all84"]; names = [str(n) for n in d["names_all84"]]; core = d["magmatic_core_mask"].astype(bool); Xbg = bg["X"].astype(float)
lab = np.load(f"{OLD}/13_Final_84x23Layers/lab_84_23_BASELINE.npy")
dep = gpd.read_file(f"{OLD}/5b_Team_Handoff/01_Deposits/Carlin_84_Training_Deposits.gpkg").to_crs(26911); dn = dep["SITE_NAME"].tolist()
pos = {n: i for i, n in enumerate(dn)}; xy84 = np.array([[dep.geometry.iloc[pos[n]].x, dep.geometry.iloc[pos[n]].y] for n in names])
FAM = {"magmatic": (lab == 1) & core, "magma_poor": lab == 0, "older_buried": lab == 2}
INTR = ["eng_eocene_belt_density", "src_eocene_intrusive_proximity", "eng_volcanic_intrusive_proximity", "eng_buried_pluton_proximity", "src_low_magmatic_fertility"]
ci = lambda ids: [IDS.index(s) for s in ids]
rng = np.random.default_rng(1); sub = rng.choice(len(Xbg), 10000, replace=False)

def cv(mask, cols, k=4, bags=8):
    X = X84[mask][:, cols]; xy = xy84[mask]; folds = KMeans(min(k, len(X) // 3), random_state=0, n_init=10).fit_predict(xy); pct = []; 
    for f in np.unique(folds):
        tr, te = folds != f, folds == f; m = fit_pu_rf(X[tr], Xbg[:, cols], n_bags=bags, n_trees=40, seed=int(f)); sr = score(m, Xbg[sub][:, cols]); st = score(m, X[te])
        pct += [float((sr < v).mean()) for v in st]
    p = np.array(pct); return round(float(p.mean()), 3), round(float((p >= 0.8).mean()), 3), len(p)

rows = []; allc = list(range(len(IDS))); noint = [j for j, s in enumerate(IDS) if s not in INTR]
for fam, mk in FAM.items():
    n = int(mk.sum()); print(fam, n, flush=True)
    for tag, cols in (("all 19 layers", allc), (f"without {len(INTR)} intrusion-dependent layers", noint)):
        mp, top20, nn = cv(mk, cols); rows.append(dict(kind="within-family spatial-block CV", train=fam, test=fam, layers=tag, n_test=nn, mean_percentile=mp, share_top20=top20))
# cross-family transfer (test family never in training)
models = {}
for fam, mk in FAM.items():
    models[fam] = fit_pu_rf(X84[mk], Xbg, n_bags=25, seed=0); joblib.dump(models[fam], f"{REF}/reference_model_{fam}.joblib", compress=3)
for a in FAM:
    sr = score(models[a], Xbg[sub])
    for b, mb in FAM.items():
        if a == b: continue
        st = score(models[a], X84[mb]); p = np.array([(sr < v).mean() for v in st]); rows.append(dict(kind="cross-family (test family unseen)", train=a, test=b, layers="all 19 layers", n_test=int(mb.sum()), mean_percentile=round(float(p.mean()), 3), share_top20=round(float((p >= 0.8).mean()), 3)))
        mdl = fit_pu_rf(X84[FAM[a]][:, noint], Xbg[:, noint], n_bags=10, seed=0); sr2 = score(mdl, Xbg[sub][:, noint]); st2 = score(mdl, X84[mb][:, noint]); p2 = np.array([(sr2 < v).mean() for v in st2])
        rows.append(dict(kind="cross-family (test family unseen)", train=a, test=b, layers=f"without {len(INTR)} intrusion-dependent layers", n_test=int(mb.sum()), mean_percentile=round(float(p2.mean()), 3), share_top20=round(float((p2 >= 0.8).mean()), 3)))
df = pd.DataFrame(rows); df.to_csv(os.path.join(ROOT, "tests", "archetype_cv_results.csv"), index=False); print(df.to_string())
np.savez_compressed(f"{REF}/reference_archetypes.npz", ids=np.array(IDS), **{f"X_{f}": X84[m] for f, m in FAM.items()}, **{f"xy_{f}": xy84[m] for f, m in FAM.items()}, **{f"names_{f}": np.array(names)[m] for f, m in FAM.items()})
meta = dict(
 magmatic=dict(n=int(FAM["magmatic"].sum()), label="Eocene-magmatic Carlin-type", use_when="Ore sits near intrusions or volcanics of the same age as mineralization (Nevada Eocene belt: dated igneous rocks within a few km of ore).", exclude_layers=[]),
 magma_poor=dict(n=int(FAM["magma_poor"].sum()), label="Magma-poor Carlin-type", use_when="Little or no magmatic input is documented: ore controlled by structure and reactive carbonate or carbonaceous host with no age-matched intrusions nearby (for example Yukon Rackla belt).", exclude_layers=INTR),
 older_buried=dict(n=int(FAM["older_buried"].sum()), label="Older-intrusion / buried-pluton Carlin-type", use_when="Mineralization is linked to older (for example Jurassic-Cretaceous) or buried intrusions rather than an age-matched belt.", exclude_layers=["eng_eocene_belt_density", "src_eocene_intrusive_proximity"]),
 caveat="Magma-poor and older/buried families have only 18 and 12 deposits. Their forests are thin; treat results as indicative and read them with the within-family cross-validation in tests/archetype_cv_results.csv.")
json.dump(meta, open(f"{REF}/archetypes.json", "w"), indent=1)
