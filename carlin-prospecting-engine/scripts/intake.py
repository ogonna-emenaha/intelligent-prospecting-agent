#!/usr/bin/env python3
"""Inventory a user's data folder and say which Carlin evidence layers can be built. READ-ONLY: never writes into the data folder.
Run this ONLY after the study-area boundary and the descriptive model have been agreed with the user (see SKILL.md).
Usage: python intake.py DATA_DIR [--json out.json]"""
import os, sys, re, json, glob, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
from layer_registry import LAYERS, status, AGREEMENT

PATHFINDERS = ["As", "Sb", "Tl", "Hg", "W", "Au", "Ag"]; IMMOBILES = ["Al", "Ti", "Th", "Zr"]; FERT = ["Sr", "Y", "La"]
# layer id -> (required roles, note). Roles are matched against what was found in the folder.
NEEDS = {
 "arch_fault_density": (["faults"], "Weighted line density, 10 km radius"),
 "arch_fault_intersection_density": (["faults"], "Needs the fault-class field: which faults are high-angle conduits and which are thrusts"),
 "arch_fluid_pathway_connectivity": (["faults"], "Total length of the connected fault network, averaged within 5 km"),
 "arch_structural_complexity": (["faults"], "1 minus mean resultant length of fault strike, 10 km radius"),
 "arch_lithologic_contact_proximity": (["geologic_map"], "Ask the user which contacts matter (for example carbonate against siliciclastic); all boundaries is a poor substitute"),
 "trap_reactive_carbonate": (["geologic_map"], "The user names the units that are reactive host rock"),
 "trap_impermeable_cap": (["geologic_map"], "The user names the units that seal"),
 "trap_argillic_k_depletion": (["radiometrics"], "Potassium lower than expected for the rock type; a lithology baseline from the geologic map is better. Another alteration proxy can be supplied as a raster"),
 "trap_fluid_chemistry": (["geochemistry"], "Pathfinder suite over immobile elements, IDW to the grid. Au/Ag version if gold assays exist"),
 "trap_preservation": (["cover_depth"], "Depth classes under cover (100/300/600 m); a cover mask from the geologic map is better"),
 "eng_eocene_belt_density": (["dated_igneous"], "The user sets the mineralization age window"),
 "eng_volcanic_intrusive_proximity": (["geologic_map"], "The user names the igneous units"),
 "eng_buried_pluton_proximity": (["gravity", "magnetics", "geologic_map"], "Gravity-low plus magnetic-high anomalies, calibrated on mapped intrusions"),
 "eng_mantle_upwelling": (["gravity"], "Gravity upward-continued 8 km"),
 "src_eocene_intrusive_proximity": (["dated_igneous"], "Same age window as the belt density layer"),
 "src_jurassic_cretaceous_intrusive_proximity": (["geologic_map"], "The user names the older intrusive units"),
 "src_basement_proximity": (["geologic_map"], "The user names the basement units; the layer is distance FROM outcrop"),
 "src_low_magmatic_fertility": (["geochemistry"], "Needs Sr, Y and La"),
 "src_crustal_gradient_fav": (["gravity"], "1 minus rescaled log total horizontal gradient of gravity"),
}

def _vector_info(path):
    import geopandas as gpd, shapely, numpy as np
    g = gpd.read_file(path); n = len(g); issue = None
    ok = g[g.geometry.notna() & ~g.geometry.is_empty]
    c = shapely.get_coordinates(ok.geometry.values) if len(ok) else np.zeros((0, 2))
    if len(c) and not np.isfinite(c).any(): issue = "ALL coordinates are infinite or NaN: the file is corrupted and cannot be used. Re-export it from the source."
    elif len(c) and (~np.isfinite(c)).any(): issue = "some geometries have non-finite coordinates; they will be dropped"
    gt = ",".join(sorted(set(ok.geom_type))) if len(ok) else "empty"
    bounds = [round(float(v)) for v in g.total_bounds] if issue is None and n else None
    return dict(kind=f"vector ({gt})", crs=(f"EPSG:{g.crs.to_epsg()}" if g.crs and g.crs.to_epsg() else str(g.crs)[:40]), n=n, columns=[x for x in g.columns if x != "geometry"], bounds=bounds, issue=issue)

def classify(path):
    import rasterio
    base = os.path.basename(path).lower(); ext = os.path.splitext(path)[1].lower(); info = dict(file=path, roles=[], issue=None, columns=[])
    try:
        if ext in (".tif", ".tiff"):
            with rasterio.open(path) as ds:
                epsg = ds.crs.to_epsg() if ds.crs else None
                info.update(kind="raster", crs=f"EPSG:{epsg}" if epsg else (ds.crs.to_string()[:30] if ds.crs else None), n=ds.count, columns=[f"{ds.width}x{ds.height}, cell {ds.res[0]:.0f}"],
                            bounds=[round(v) for v in ds.bounds])
            if re.search(r"grav|bouguer", base): info["roles"].append("gravity")
            elif re.search(r"mag|aeromag|tmi|rtp", base) and "magma" not in base: info["roles"].append("magnetics")
            elif re.search(r"(^|[_\- ])k([_\-. ]|$)|potassium|radiom", base): info["roles"].append("radiometrics")
            elif re.search(r"(^|[_\- ])th([_\-. ]|$)|thorium", base): info["roles"].append("radiometrics_th")
            elif re.search(r"depth|basement|cover|overburden", base): info["roles"].append("cover_depth")
            elif re.search(r"dem|elev", base): info["roles"].append("dem")
            if info["crs"] and info["crs"] != "EPSG:26911" and "NAD_1983_Albers" in str(info["crs"]): info["note"] = "different projection from the other layers; build_layers reprojects it"
        elif ext in (".gpkg", ".shp", ".geojson", ".json"):
            info.update(_vector_info(path)); cols = " ".join(info["columns"]).lower(); kind = info["kind"].lower()
            if "line" in kind and re.search(r"fault|lineament|structur|thrust|shear", base + " " + cols): info["roles"].append("faults")
            if "polygon" in kind and re.search(r"geol|lith|bedrock|map.?unit|formation", base + " " + cols): info["roles"].append("geologic_map")
            if "point" in kind and re.search(r"\bas_ppm|\bsb_ppm|\btl_ppm|geochem|assay|stream|sediment", base + " " + cols): info["roles"].append("geochemistry")
            if "point" in kind and re.search(r"deposit|mine|occurr|prospect|training", base): info["roles"].append("deposits")
            if "point" in kind and re.search(r"\bage\b|age_ma", cols): info["roles"].append("dated_igneous")
            if re.search(r"boundar|tract|study.?area|aoi", base): info["roles"].append("boundary")
        elif ext == ".csv":
            import pandas as pd
            d = pd.read_csv(path, nrows=50, low_memory=False); info.update(kind="table", n=None, columns=list(d.columns)[:30], crs=None)
            if re.search(r"ar40ar39|kar_|samples_usgs|geochron", base): info["roles"].append("geochron_table")
        else: return None
    except Exception as e:
        info["issue"] = f"could not be read: {e}"; info.setdefault("kind", "unreadable")
    return info

def scan(folder):
    files = sorted({f for ext in ("*.tif", "*.tiff", "*.gpkg", "*.shp", "*.geojson", "*.csv") for f in glob.glob(os.path.join(folder, "**", ext), recursive=True)})
    return [c for c in (classify(f) for f in files) if c]

def roles_found(inv):
    roles = {r for i in inv if not i.get("issue") or "dropped" in (i.get("issue") or "") for r in i["roles"]}
    gc = [i for i in inv if "geochemistry" in i["roles"] and not i.get("issue")]
    names = {i["file"] for i in inv if "geochron_table" in i["roles"]}
    if sum(any(k in os.path.basename(f).lower() for k in ks) for ks in (("ar40ar39",), ("kar_",), ("samples_usgs",)) for f in names) >= 3 or \
       (any("samples" in os.path.basename(f).lower() for f in names) and len(names) >= 2): roles.add("dated_igneous")
    if "radiometrics_th" in roles: pass
    return roles, gc

def geochem_elements(gc):
    out = {}
    for i in gc:
        cols = {c.split("_")[0].lower(): c for c in i["columns"]}
        out[i["file"]] = dict(pathfinders=[e for e in PATHFINDERS if e.lower() in cols], immobiles=[e for e in IMMOBILES if e.lower() in cols], fertility=[e for e in FERT if e.lower() in cols])
    return out

def plan(inv):
    roles, gc = roles_found(inv); ge = geochem_elements(gc); out = []
    best = {"pathfinders": set(), "immobiles": set(), "fertility": set()}
    for v in ge.values():
        for k in best: best[k] |= set(v[k])
    for lid, name, part, cat, inputs in LAYERS:
        need, note = NEEDS[lid]; missing = [r for r in need if r not in roles]
        extra = ""
        if lid == "trap_fluid_chemistry" and not missing and not (len(best["pathfinders"] & {"As", "Sb", "Tl", "Hg"}) >= 2 and len(best["immobiles"]) >= 2): missing = ["geochemistry with at least two pathfinders and two immobile elements"]
        if lid == "src_low_magmatic_fertility" and not missing and best["fertility"] != set(FERT): missing = ["Sr, Y and La in the geochemistry"]
        out.append(dict(layer=name, id=lid, part=part, state="can build" if not missing else "missing inputs", missing=missing, builder_status=status(lid), agreement=AGREEMENT.get(lid), note=note))
    return roles, out, ge

if __name__ == "__main__":
    folder = sys.argv[1]; inv = scan(folder); roles, pl, ge = plan(inv)
    print("(Reminder: confirm the study-area boundary and the descriptive model with the user before building anything.)\n")
    print(f"Files found: {len(inv)}\n")
    for i in inv:
        print(f"- {os.path.relpath(i['file'], folder)} | {i.get('kind')} | {i.get('crs')} | looks like: {', '.join(i['roles']) or 'unclassified (ask the user)'}")
        if i.get("issue"): print(f"    PROBLEM: {i['issue']}")
        if i.get("note"): print(f"    note: {i['note']}")
    if ge:
        for f, v in ge.items(): print(f"\nGeochemistry elements found in {os.path.basename(f)}: pathfinders {v['pathfinders']}, immobiles {v['immobiles']}, fertility {v['fertility']}")
    nb = sum(p["state"] == "can build" for p in pl); print(f"\nCan build {nb} of {len(pl)} layers with the detected data. 'substitute' means the builder does not reproduce the Nevada layer, so the user must accept it as a proxy.\n")
    for p in pl:
        miss = f" (needs {', '.join(p['missing'])})" if p["missing"] else ""
        print(f"[{p['state']:14s}] {p['part']:12s} {p['layer']} | {p['builder_status']}{miss}")
    if "--json" in sys.argv: json.dump(dict(inventory=inv, roles=sorted(roles), plan=pl), open(sys.argv[sys.argv.index("--json") + 1], "w"), indent=1, default=str)
