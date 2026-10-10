#!/usr/bin/env python3
"""Build the user's evidence layers from a JSON config that the user has approved layer by layer. Never writes into the input data folders.
Usage: python build_layers.py config.json OUT_DIR
A study-area boundary is REQUIRED (config 'tract' with a polygon file, or explicit bounds). Every output is a GeoTIFF named <layer_id>.tif,
0-1, nodata -9999, rescaled within the study area (Nevada convention). See references/example_config.json."""
import os, sys, json, re, warnings; warnings.filterwarnings("ignore")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import builders as B
from io_utils import write_raster
from layer_registry import status, AGREEMENT

_CACHE = {}
def _vec(path, crs):
    import geopandas as gpd
    if (path, crs) not in _CACHE:
        g = gpd.read_file(path); g = g.to_crs(crs) if g.crs else g.set_crs(crs); _CACHE[(path, crs)] = B.clean_geoms(g)
    return _CACHE[(path, crs)]

PATH_KEYS = ("file", "raster", "gravity", "magnetics", "k", "th", "depth", "ar", "kar", "samples")
def _resolve(obj, base):
    """Join relative paths in the config to cfg['data_dir'] (read-only data folder)."""
    if isinstance(obj, dict): return {k: (os.path.join(base, v) if k in PATH_KEYS and isinstance(v, str) and not os.path.isabs(v) else _resolve(v, base)) for k, v in obj.items()}
    return obj

def _match(series, regex): return series.fillna("").astype(str).str.contains(regex, case=False, regex=True)

def _warp(path, grid, crs):
    import rasterio, rasterio.warp
    with rasterio.open(path) as src:
        a = np.full((grid.H, grid.W), np.nan)
        rasterio.warp.reproject(rasterio.band(src, 1), a, src_transform=src.transform, src_crs=src.crs, dst_transform=grid.transform, dst_crs=crs, resampling=rasterio.warp.Resampling.bilinear, dst_nodata=np.nan, src_nodata=src.nodata)
    return np.where(a < -9000, np.nan, a)

def _unit_raster(grid, g, field):
    import rasterio.features
    f = g[field].fillna("").astype(str); codes = {c: i + 1 for i, c in enumerate(sorted(f.unique()))}
    return rasterio.features.rasterize([(x, codes[c]) for x, c in zip(g.geometry, f) if x is not None and not x.is_empty], out_shape=(grid.H, grid.W), transform=grid.transform, fill=0, dtype="int32")

def _dated(c, crs):
    if "geochron" in c:
        gc = c["geochron"]; return B.load_geochron(gc.get("ar"), gc.get("kar"), gc["samples"], crs=crs)
    g = _vec(c["file"], crs); g = g.rename(columns={c["age_field"]: "age"}); g["age"] = g["age"].astype(float); return g

SCALE_FILE = os.path.join(os.path.dirname(__file__), "..", "reference", "nevada_scale.json")
def _load_scale(mode):
    """mode 'fixed' (default): cut-offs calibrated on the Nevada tract, identical in every region. 'relative': legacy, rescaled inside the study area."""
    if mode == "relative": return {}
    if not os.path.exists(SCALE_FILE): raise FileNotFoundError("reference/nevada_scale.json is missing; run tools/calibrate_nevada_scale.py or set scaling to 'relative'")
    return json.load(open(SCALE_FILE))

def _harmonise_faults(g, c):
    """Optional mapping-scale control for fault layers: config keys scale_field (map-scale denominator in thousands) and common_scale (default 250)."""
    if not c.get("scale_field"): return list(g.geometry), None
    geoms, keep = B.generalise_lines(list(g.geometry), g[c["scale_field"]].values, c.get("common_scale", 250), c.get("simplify_m", 125.0), c.get("min_length_m", 2000.0))
    gg = [x for x, k in zip(geoms, keep) if k]
    return gg, dict(features_in=len(geoms), features_kept=int(keep.sum()), common_scale=c.get("common_scale", 250))

def build(cfg, out_dir):
    import rasterio
    if "tract" not in cfg or not cfg["tract"]: raise ValueError("No study-area boundary. Ask the user for the boundary (polygon file or bounds) before building.")
    mode = cfg.get("scaling", "fixed"); SC = _load_scale(mode)
    cfg = _resolve(cfg, cfg.get("data_dir", "")) if cfg.get("data_dir") else cfg
    os.makedirs(out_dir, exist_ok=True); crs = cfg["crs"]; cell = float(cfg.get("cell", 1000))
    if "file" in cfg["tract"]:
        tg = _vec(cfg["tract"]["file"], crs); bounds = tuple(cfg["tract"].get("bounds", tg.total_bounds)); tract_geom = list(tg.geometry)
    elif "bounds" in cfg["tract"]: bounds = tuple(cfg["tract"]["bounds"]); tract_geom = None
    else: raise ValueError("tract needs 'file' or 'bounds'")
    grid = B.Grid(bounds, cell, crs); prof = dict(transform=grid.transform, crs=rasterio.crs.CRS.from_string(crs), shape=(grid.H, grid.W))
    tract = B.polygons_to_mask(grid, tract_geom) if tract_geom else np.ones((grid.H, grid.W), bool)
    done = {}
    for lid, c in cfg["layers"].items():
        print("building", lid, flush=True)
        B.LAST.clear(); fx = SC.get(lid, {}); fixed = tuple(fx["fixed"]) if "fixed" in fx else None; cap = fx.get("cap_m"); harm = None
        try:
            if "raster" in c:
                a = _warp(c["raster"], grid, crs); arr = B.minmax(a, tract & np.isfinite(a)); arr = 1 - arr if c.get("invert") else arr
            elif lid in ("arch_fault_density", "arch_fault_intersection_density", "arch_structural_complexity", "arch_fluid_pathway_connectivity"):
                g = _vec(c["file"], crs)
                if c.get("scale_field") and not str(c["scale_field"]).startswith("(optional)"):
                    # mapping-scale control: put every fault on one common mapping detail before any fault layer is computed
                    geoms_all, keep = B.generalise_lines(list(g.geometry), g[c["scale_field"]].values, c.get("common_scale", 250), c.get("simplify_m", 125.0), c.get("min_length_m", 2000.0))
                    g = g.assign(geometry=geoms_all)[keep]; harm = dict(features_in=int(len(keep)), features_kept=int(keep.sum()), common_scale=c.get("common_scale", 250))
                geoms = list(g.geometry)
                if lid == "arch_fault_density":
                    w = None
                    if c.get("weight_field") and c.get("weights"): w = g[c["weight_field"]].map(lambda t: next((v for k, v in c["weights"].items() if k.lower() in str(t).lower()), 1.0)).values
                    arr = B.fault_density(grid, tract, geoms, w, c.get("radius_m", 10000), fixed)
                elif lid == "arch_fault_intersection_density":
                    if c.get("field") and c.get("a_regex") and c.get("b_regex"):
                        A = list(g[_match(g[c["field"]], c["a_regex"])].geometry); Bq = list(g[_match(g[c["field"]], c["b_regex"])].geometry); arr = B.fault_intersection_proximity(grid, tract, A, Bq, cap)
                    else: arr = B.fault_intersection_proximity(grid, tract, geoms, None, cap)
                elif lid == "arch_structural_complexity": arr = B.structural_complexity(grid, tract, geoms, c.get("radius_m", 10000), fixed=fixed)
                else: arr = B.fault_network_connectivity(grid, tract, geoms, c.get("radius_km", 5), fixed)
            elif lid in ("eng_eocene_belt_density", "src_eocene_intrusive_proximity"):
                g = _dated(c, crs); x, y, age = g.geometry.x.values, g.geometry.y.values, g["age"].values
                arr = B.age_window_density(grid, tract, x, y, age, c["age_min"], c["age_max"], c.get("radius_m", 140000), fixed) if lid == "eng_eocene_belt_density" else B.age_window_proximity(grid, tract, x, y, age, c["age_min"], c["age_max"], cap)
            elif lid in ("eng_volcanic_intrusive_proximity", "trap_reactive_carbonate", "trap_impermeable_cap", "src_jurassic_cretaceous_intrusive_proximity", "src_basement_proximity", "arch_lithologic_contact_proximity"):
                g = _vec(c["file"], crs); sel = _match(g[c["field"]], c["regex"])
                if c.get("age_field"): sel &= _match(g[c["age_field"]], c["age_regex"])
                polys = list(g[sel].geometry)
                if not polys: raise ValueError("no geologic units matched the text given; ask the user to confirm the unit names")
                if lid == "src_basement_proximity": arr = B.basement_distance(grid, tract, polys, cap)
                elif lid == "arch_lithologic_contact_proximity":
                    is_line = bool(g[sel].geom_type.str.contains("Line").all())   # contacts supplied as lines, or unit polygons whose boundaries are the contacts
                    arr = B.proximity_score(B.lines_to_mask(grid, polys) if is_line else B.boundary_mask(grid, polys), tract, cap_m=cap)
                else: arr = B.proximity_score(B.polygons_to_mask(grid, polys), tract, cap_m=cap)
            elif lid in ("eng_mantle_upwelling", "src_crustal_gradient_fav"):
                a = _warp(c["gravity"], grid, crs); arr = B.upward_continuation(a, tract, c.get("height_m", 8000), cell) if lid == "eng_mantle_upwelling" else B.total_horizontal_gradient(a, tract, cell)
            elif lid == "eng_buried_pluton_proximity":
                ig = c["intrusive"]; g = _vec(ig["file"], crs); im = B.polygons_to_mask(grid, list(g[_match(g[ig["field"]], ig["regex"])].geometry)) & tract
                arr = B.buried_pluton_proximity(grid, tract, _warp(c["gravity"], grid, crs), _warp(c["magnetics"], grid, crs), im)
            elif lid == "trap_argillic_k_depletion":
                K = _warp(c["k"], grid, crs)
                if c.get("lithology"): lg = _vec(c["lithology"]["file"], crs); arr = B.potassium_depletion(grid, tract, K, _unit_raster(grid, lg, c["lithology"]["field"]))
                else: arr = B.potassium_depletion_th(grid, tract, K, _warp(c["th"], grid, crs))
            elif lid == "trap_fluid_chemistry":
                g = _vec(c["file"], crs)
                arr = B.gold_silver_footprint(grid, tract, g, c.get("au", "Au_ppm"), c.get("ag", "Ag_ppm")) if c.get("mode") == "au_ag" else B.geochem_pathfinder_score(grid, tract, g, c.get("pathfinders", ["As_ppm", "Sb_ppm", "Tl_ppm"]), c.get("immobiles", ["Al_pct", "Ti_pct", "Th_ppm"]))
            elif lid == "src_low_magmatic_fertility":
                arr = B.magmatic_fertility(grid, tract, _vec(c["file"], crs), c.get("sr", "Sr_ppm"), c.get("y", "Y_ppm"), c.get("la", "La_ppm"))
            elif lid == "trap_preservation":
                D = _warp(c["depth"], grid, crs); cm = None
                if c.get("cover"): cg = _vec(c["cover"]["file"], crs); cm = B.polygons_to_mask(grid, list(cg[_match(cg[c["cover"]["field"]], c["cover"]["regex"])].geometry))
                arr = B.preservation_from_depth(grid, tract, D, cm, tuple(c.get("breaks", (100, 300, 600))))
            else: raise NotImplementedError(f"no builder for {lid}; supply a prepared raster with {{'raster': path}} or leave the layer out")
            if arr is None: raise ValueError("builder returned no data")
            arr = np.where(tract, arr, np.nan); write_raster(os.path.join(out_dir, lid + ".tif"), arr, prof)
            v = arr[tract & np.isfinite(arr)]; fixed_scale = bool(fx) and "raster" not in c
            done[lid] = dict(result="built", proxy_status=("user-supplied raster (equivalence not tested)" if "raster" in c else status(lid)), agreement_with_nevada_layer=None if "raster" in c else AGREEMENT.get(lid),
                             scaling=("fixed (Nevada-calibrated physical cut-offs)" if fixed_scale else "area-relative (rescaled inside this study area; values are not comparable across regions)"),
                             cutoffs_used={k: float(f"{float(val):.9g}") for k, val in B.LAST.items()},
                             share_at_zero=round(float((v <= 1e-9).mean()), 3), share_at_one=round(float((v >= 1 - 1e-9).mean()), 3), mapping_scale_control=harm)
        except Exception as e:
            done[lid] = dict(result=f"FAILED: {e}")
    meta = dict(scaling_mode=mode, cell_m=cell, study_area_km2=round(float(tract.sum()) * cell * cell / 1e6, 1),
                tract_declared_permissive=bool(cfg.get("tract_is_permissive", False)), tract_basis=cfg.get("tract_basis"),
                warning=None if cfg.get("tract_is_permissive", False) else "The study area was not declared a permissive tract (ground where Carlin-type deposits can form). Cells outside permissive ground make the comparison unlike-for-like.")
    done["_meta"] = meta
    json.dump(done, open(os.path.join(out_dir, "build_log.json"), "w"), indent=1); return done

if __name__ == "__main__":
    cfg = json.load(open(sys.argv[1])); print(json.dumps(build(cfg, sys.argv[2]), indent=1))
