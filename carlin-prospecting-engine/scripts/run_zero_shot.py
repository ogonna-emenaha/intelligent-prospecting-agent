#!/usr/bin/env python3
"""Zero-shot Carlin transfer: apply a frozen Nevada archetype to a new district. No local deposits needed.
Usage: python run_zero_shot.py --layers DIR --out DIR --archetype magmatic|magma_poor|older_buried [--deposits points.gpkg|csv for VALIDATION ONLY]
DIR holds GeoTIFFs named <layer_id>.tif (see layer_registry.py), all on one grid, values 0-1, nodata -9999 or nan.
Requires a permissive tract (build_log.json written by build_layers.py records it; otherwise pass --tract-basis TEXT or --allow-unfiltered)."""
import os, sys, json, argparse, warnings; warnings.filterwarnings("ignore")
import numpy as np, joblib
sys.path.insert(0, os.path.dirname(__file__))
from engine import fit_pu_rf, score, domain_model, in_domain, percentile_rank, layer_shift_table, archetype_fidelity
from layer_registry import IDS, PART, NAME, status, AGREEMENT, FIXED_SCALE, DIAGNOSTIC
from io_utils import read_layers, write_raster, deposits_xy
REF = os.path.join(os.path.dirname(__file__), "..", "reference")
PARTS = ["Source", "Energy", "Pathway", "Trap", "Preservation"]

def load_archetype(name):
    meta = json.load(open(os.path.join(REF, "archetypes.json"))); arch = np.load(os.path.join(REF, "reference_archetypes.npz"), allow_pickle=True)
    if name not in meta or name == "caveat": raise ValueError(f"unknown archetype {name!r}; choose one of {[k for k in meta if k != 'caveat']}")
    return meta[name], arch[f"X_{name}"], arch[f"xy_{name}"], [str(x) for x in arch["ids"]], meta["caveat"]

def run(layers_dir, out_dir, min_layers=10, min_parts=3, top_share=0.20, seed=0, n_bags=25, deposits_xy=None, exclude=(), archetype="magmatic",
        tract_basis=None, allow_unfiltered=False, diagnostic_weight=0.0, keep_archetype_exclusions=False):
    os.makedirs(out_dir, exist_ok=True)
    bl = os.path.join(layers_dir, "build_log.json"); blog = json.load(open(bl)) if os.path.exists(bl) else {}; meta_b = blog.get("_meta", {})
    declared = bool(meta_b.get("tract_declared_permissive")) or bool(tract_basis)
    am, Xarch, xyarch, aids, caveat = load_archetype(archetype)
    rep = dict(archetype=dict(name=archetype, label=am["label"], n_deposits=int(len(Xarch)), use_when=am["use_when"], caveat=caveat),
               comparison_basis=dict(tract_declared_permissive=declared, tract_basis=tract_basis or meta_b.get("tract_basis"), study_area_km2=meta_b.get("study_area_km2")))
    if not declared and not allow_unfiltered:
        rep["verdict"] = "Not run: permissive tract required"
        rep["reason"] = ("A zero-shot run compares the study area with Nevada ground that sits inside a permissive tract (rock and structure where Carlin-type deposits can form). "
                         "Running it over a whole territory that includes terranes where Carlin deposits do not form is an unlike-for-like comparison. Supply the permissive tract as the study boundary "
                         "(build_layers: tract_is_permissive=true and tract_basis=<why this ground is permissive>), pass --tract-basis, or pass --allow-unfiltered to run anyway and have it recorded.")
        json.dump(rep, open(os.path.join(out_dir, "zero_shot_summary.json"), "w"), indent=1); return rep
    if not declared: rep["comparison_basis"]["warning"] = "RUN WITHOUT A PERMISSIVE TRACT (--allow-unfiltered). Cells on non-permissive ground are scored against Nevada permissive ground; treat the map as unlike-for-like."
    excl_arch = [] if keep_archetype_exclusions else [s for s in am["exclude_layers"] if s not in exclude]
    L, prof = read_layers(layers_dir, IDS); skip = set(exclude) | set(excl_arch); have = [s for s in IDS if s in L and s not in skip]
    parts_have = sorted({PART[s] for s in have}); missing = [s for s in IDS if s not in L]
    fixed_flags = {s: ("fixed" in blog.get(s, {}).get("scaling", "")) if blog else (s in FIXED_SCALE) for s in have}
    rep.update(layers_excluded_by_user=[NAME[s] for s in exclude if s in NAME], layers_excluded_by_archetype=[NAME[s] for s in excl_arch if s in L], layers_available=[NAME[s] for s in have],
               proxy_layers=[dict(layer=NAME[s], status=status(s), agreement_with_nevada_layer=AGREEMENT.get(s)) for s in have if status(s) != "matches the Nevada layer"], layers_missing=[NAME[s] for s in missing],
               parts_represented=parts_have, parts_missing=[p for p in PARTS if p not in parts_have],
               layers_area_relative=[NAME[s] for s in have if not fixed_flags[s]],
               scaling_note="Layers marked area-relative are rescaled inside the study area, so a value of 0.8 does not mean the same thing here as in Nevada. Fixed-scale layers use Nevada-calibrated physical cut-offs (km, km per km2) and are comparable.")
    diag = [s for s in have if s in DIAGNOSTIC]; ctx = [s for s in have if s not in DIAGNOSTIC]
    rep["diagnostic_evidence"] = dict(layers=[NAME[s] for s in diag], present=bool(diag),
        note=None if diag else "No ore-fluid chemistry or alteration layer is present. The map shows where host rock, structure and intrusions resemble Nevada Carlin-type ground; those exist in most sedimentary terranes, so it cannot tell Carlin-type from other sediment-hosted systems. Add pathfinder geochemistry (As, Sb, Tl, Hg, Au) or an alteration layer to make it Carlin-specific.")
    core_run = len(have) >= 8 and "Pathway" in parts_have and "Trap" in parts_have   # transferable core: Pathway + Trap (+ diagnostic chemistry)
    if core_run and not (len(have) >= min_layers and len(parts_have) >= min_parts):
        rep["coverage_note"] = ("Core run: Pathway and Trap are represented; " + ", ".join(p for p in ("Source", "Energy", "Preservation") if p not in parts_have) +
                                " not represented. Pathway and Trap are the parts that transfer from Nevada; Source and Energy are district-specific.")
    if not core_run and (len(have) < min_layers or len(parts_have) < min_parts):
        rep["verdict"] = "Insufficient evidence"
        rep["reason"] = f"{len(have)} layers across {len(parts_have)} of 5 mineral-system parts; need at least {min_layers} layers and {min_parts} parts (provisional defaults)."
        json.dump(rep, open(os.path.join(out_dir, "zero_shot_summary.json"), "w"), indent=1); return rep
    bg = np.load(os.path.join(REF, "reference_background.npz"), allow_pickle=True)
    rid = [str(x) for x in bg["ids"]]; cols = [rid.index(s) for s in have]; acols = [aids.index(s) for s in have]
    Xdep = Xarch[:, acols]; Xbg = bg["X"][:, cols].astype(float)
    tract = np.all(np.isfinite(np.stack([L[s] for s in have])), 0)
    X = np.stack([L[s][tract] for s in have], 1)
    model = fit_pu_rf(Xdep, Xbg, n_bags=n_bags, seed=seed); prob = score(model, X); pr = percentile_rank(prob)
    if diagnostic_weight > 0 and diag and ctx:
        ci = [have.index(s) for s in ctx]; di = [have.index(s) for s in diag]
        m_ctx = fit_pu_rf(Xdep[:, ci], Xbg[:, ci], n_bags=n_bags, seed=seed); m_dia = fit_pu_rf(Xdep[:, di], Xbg[:, di], n_bags=n_bags, seed=seed)
        pr = (1 - diagnostic_weight) * percentile_rank(score(m_ctx, X[:, ci])) + diagnostic_weight * percentile_rank(score(m_dia, X[:, di])); pr = percentile_rank(pr); prob = pr
        rep["diagnostic_evidence"]["weight_in_map"] = diagnostic_weight
    dm = domain_model(Xbg); ok, dist = in_domain(dm, X)
    def grid(v): g = np.full(tract.shape, np.nan); g[tract] = v; return g
    write_raster(os.path.join(out_dir, "carlin_similarity_probability.tif"), grid(prob), prof)
    write_raster(os.path.join(out_dir, "carlin_similarity_percentile.tif"), grid(pr), prof)
    cls = np.where(~ok, 0, np.where(pr >= 1 - top_share, 2, 1))   # 0 outside Nevada reference domain, 1 in domain, 2 in domain and top share
    write_raster(os.path.join(out_dir, "reference_domain_class.tif"), grid(cls.astype(float)), prof)
    out_share = float((~ok).mean()); shifts = layer_shift_table(Xbg, X, have)
    rep.update(n_cells=int(tract.sum()), share_outside_reference_domain=round(out_share, 3), layers_to_check=[r for r in shifts if r["flag"]], layer_shift_table=shifts, share_in_domain_top=round(float((cls == 2).mean()), 3),
               verdict=("Outside reference domain" if out_share > 0.5 else "Carlin-like evidence pattern present (not a deposit probability)"),
               notes=["The reference-domain flag is experimental: it catches gross extrapolation only and has not been shown to separate Carlin from non-Carlin ground.",
                      f"Archetype '{archetype}' = {len(Xarch)} Nevada deposits; the forest is fitted on the available layers only.",
                      "Provisional coverage rule: " + f"{min_layers} layers and {min_parts} parts."])
    if deposits_xy is not None and len(deposits_xy):
        import rasterio
        r, c = rasterio.transform.rowcol(prof["transform"], deposits_xy[:, 0], deposits_xy[:, 1])
        r = np.array(r); c = np.array(c); okp = (r >= 0) & (r < tract.shape[0]) & (c >= 0) & (c < tract.shape[1])
        pg = grid(pr); vals = pg[r[okp], c[okp]]; vals = vals[np.isfinite(vals)]
        rep["validation"] = dict(n_deposits_in_area=int(len(vals)), mean_percentile=round(float(vals.mean()), 3) if len(vals) else None,
                                 share_in_top_20pct=round(float((vals >= 0.8).mean()), 3) if len(vals) else None, chance_share=0.2,
                                 note="Quick check only. For a proper test use run_evaluation.py (area hold-out, baseline, negative control).")
    try:
        full = joblib.load(os.path.join(REF, f"reference_model_{archetype}.joblib"))
        rep["archetype_comparison"] = archetype_fidelity(Xarch, bg["X"].astype(float), xyarch, [aids.index(s) for s in have], full) if len(have) < len(IDS) else dict(all_layers_present=True)
    except Exception as e:
        rep["archetype_comparison"] = dict(error=str(e))
    json.dump(rep, open(os.path.join(out_dir, "zero_shot_summary.json"), "w"), indent=1); return rep

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--layers", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--archetype", default="magmatic", choices=["magmatic", "magma_poor", "older_buried"]); ap.add_argument("--deposits")
    ap.add_argument("--min-layers", type=int, default=10); ap.add_argument("--min-parts", type=int, default=3); ap.add_argument("--exclude", default="", help="comma-separated layer ids to leave out")
    ap.add_argument("--tract-basis", help="why the study area is permissive ground (use when layers were not built with build_layers)"); ap.add_argument("--allow-unfiltered", action="store_true")
    ap.add_argument("--diagnostic-weight", type=float, default=0.0, help="0-1 share of the map taken from diagnostic evidence (ore-fluid chemistry, alteration); 0 = single joint forest")
    ap.add_argument("--keep-archetype-exclusions", action="store_true", help="keep the layers the archetype normally drops (intrusion-dependent layers for magma-poor)")
    a = ap.parse_args(); xy = deposits_xy(a.deposits, a.layers) if a.deposits else None
    print(json.dumps(run(a.layers, a.out, a.min_layers, a.min_parts, deposits_xy=xy, exclude=[x for x in a.exclude.split(',') if x], archetype=a.archetype, tract_basis=a.tract_basis,
                         allow_unfiltered=a.allow_unfiltered, diagnostic_weight=a.diagnostic_weight, keep_archetype_exclusions=a.keep_archetype_exclusions), indent=1))
