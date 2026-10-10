"""Core Carlin engine: frozen Nevada reference (magmatic family), zero-shot scoring, local training.
Recipe follows the Terra-X Nevada model: standardise -> PCA (90% variance) -> PU-RandomForest (25 bags x 60 trees)."""
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier

def fit_pu_rf(Xdep, Xbg, n_bags=25, n_trees=60, seed=0, pca_fit_set=None, var=0.90):
    """Positive-unlabeled random forest. Xdep: positives (n x k). Xbg: unlabeled background sample (m x k).
    Default: forest on the layers directly (no PCA). Spatial-block tests showed the Nevada recipe's PCA is
    unstable when fitted on one family's deposits only (map correlation 0.25 vs 0.80 when fitted on all 84),
    so PCA is used only when pca_fit_set (a multi-family or background matrix) is supplied."""
    if pca_fit_set is None:
        sc = pca = None; tf = lambda X: X
    else:
        sc = StandardScaler().fit(pca_fit_set); pca = PCA(var, random_state=0).fit(sc.transform(pca_fit_set))
        tf = lambda X: pca.transform(sc.transform(X))
    P = tf(Xdep); B = tf(Xbg)
    rng = np.random.default_rng(seed); forests = []
    for _ in range(n_bags):
        ni = rng.choice(len(B), len(P), replace=False)
        X = np.vstack([P, B[ni]]); y = np.r_[np.ones(len(P)), np.zeros(len(P))]
        forests.append(RandomForestClassifier(n_trees, n_jobs=-1, random_state=int(rng.integers(1e6))).fit(X, y))
    return {"scaler": sc, "pca": pca, "forests": forests, "n_dep": len(P), "k": Xdep.shape[1]}

def _tf(model, X):
    return X if model["pca"] is None else model["pca"].transform(model["scaler"].transform(X))

def score(model, X):
    Z = _tf(model, X)
    return np.mean([f.predict_proba(Z)[:, 1] for f in model["forests"]], axis=0)

def percentile_rank(v):
    return np.argsort(np.argsort(v)) / (len(v) - 1)

def domain_model(Xref_bg, n_ref=20000, q=99.9, seed=0):
    """Extrapolation flag: robust Mahalanobis distance to the Nevada reference background (MinCovDet, 90% support).
    Experimental. Tests on Nevada layers rebuilt from raw data showed a kNN distance flag was far too sensitive to small
    builder differences (76% of Nevada flagged, same as column-shuffled data), so a robust Mahalanobis flag is used instead.
    It catches gross extrapolation only; it has NOT been shown to separate Carlin from non-Carlin ground."""
    from sklearn.covariance import MinCovDet
    rng = np.random.default_rng(seed); idx = rng.choice(len(Xref_bg), min(n_ref, len(Xref_bg)), replace=False)
    R = Xref_bg[idx]; sc = StandardScaler().fit(R); Rz = sc.transform(R)
    mcd = MinCovDet(support_fraction=0.9, random_state=0).fit(Rz)
    return {"scaler": sc, "mcd": mcd, "threshold": float(np.percentile(mcd.mahalanobis(Rz), q)), "q": q}

def in_domain(dm, X):
    d = dm["mcd"].mahalanobis(dm["scaler"].transform(X))
    return d <= dm["threshold"], d

def layer_shift_table(Xref_bg, X, ids):
    """Per-layer check: where does the new district's median fall inside the Nevada background distribution?
    Flags layers whose median sits outside Nevada's 5th-95th percentile range (check how the layer was built)."""
    rows = []
    for j, s in enumerate(ids):
        lo, hi = np.percentile(Xref_bg[:, j], [5, 95]); med = float(np.median(X[:, j]))
        rows.append(dict(layer=s, district_median=round(med, 3), nevada_p5=round(float(lo), 3), nevada_p50=round(float(np.median(Xref_bg[:, j])), 3), nevada_p95=round(float(hi), 3), flag=bool(med < lo or med > hi)))
    return rows

def signature_table(Xdep, Xbg, ids):
    """Geological signature: deposit vs background value distribution per layer."""
    import pandas as pd
    from scipy.stats import mannwhitneyu
    from sklearn.metrics import roc_auc_score
    rows = []
    for j, s in enumerate(ids):
        a, b = Xdep[:, j], Xbg[:, j]
        p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        auc = roc_auc_score(np.r_[np.ones(len(a)), np.zeros(len(b))], np.r_[a, b])
        rows.append(dict(layer=s, deposit_mean=a.mean(), background_mean=b.mean(), delta=a.mean() - b.mean(),
                         dep_p05=np.percentile(a, 5), dep_p50=np.median(a), dep_p95=np.percentile(a, 95),
                         AUC=auc, p_value=p))
    return pd.DataFrame(rows)

def archetype_fidelity(ref_dep, ref_bg, ref_xy, cols, full_model, n_bags=6, n_trees=40, seed=0):
    """DIAGNOSTIC ONLY. Compares a model built on the user's layer subset (`cols`) with the full Carlin-type archetype, both applied to Nevada.
    It does not change the user's map and does not measure its accuracy. Returns agreement of the two Nevada maps and the held-out
    Nevada deposit ranking of each model under spatial-block cross-validation."""
    from scipy.stats import spearmanr
    from sklearn.cluster import KMeans
    rng = np.random.default_rng(seed); k_all = ref_dep.shape[1]
    sub = rng.choice(len(ref_bg), min(10000, len(ref_bg)), replace=False); Xs = ref_bg[sub]
    if len(cols) == k_all: return dict(all_layers_present=True)
    reduced = fit_pu_rf(ref_dep[:, cols], ref_bg[:, cols], n_bags=25, seed=seed)
    s_full = score(full_model, Xs); s_red = score(reduced, Xs[:, cols])
    rho = float(spearmanr(s_full, s_red).correlation)
    top = lambda v: v >= np.quantile(v, 0.90); overlap = float((top(s_full) & top(s_red)).sum() / max(top(s_full).sum(), 1))
    k = min(5, len(ref_dep) // 4); folds = KMeans(k, random_state=0, n_init=10).fit_predict(ref_xy); res = {"full": [], "reduced": []}
    for f in range(k):
        tr, te = folds != f, folds == f
        for name, cc in (("full", list(range(k_all))), ("reduced", list(cols))):
            m = fit_pu_rf(ref_dep[tr][:, cc], ref_bg[:, cc], n_bags=n_bags, n_trees=n_trees, seed=f)
            sr = score(m, Xs[:, cc]); st = score(m, ref_dep[te][:, cc]); res[name] += [float((sr < v).mean()) for v in st]
    summ = lambda p: dict(auc=round(float(np.mean(p)), 3), share_top5=round(float((np.array(p) >= 0.95).mean()), 3))
    return dict(all_layers_present=False, layers_used=len(cols), layers_in_archetype=k_all, nevada_map_agreement_spearman=round(rho, 3),
                top10_overlap_with_full_archetype=round(overlap, 3), nevada_heldout_deposits_full_archetype=summ(res["full"]), nevada_heldout_deposits_your_layers=summ(res["reduced"]),
                note=("Diagnostic only. This comparison does not change your map and is not a measure of its accuracy. It shows how closely a model built on your "
                      f"{len(cols)} layers matches the full {k_all}-layer Carlin-type archetype when both are applied to Nevada. Lower agreement means the layers you "
                      "do not have carry information about the archetype that this run cannot use."))
