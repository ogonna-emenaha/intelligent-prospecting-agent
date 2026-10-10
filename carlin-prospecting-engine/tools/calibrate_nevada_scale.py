#!/usr/bin/env python3
"""Calibrate the fixed physical cut-offs used by build_layers (scaling='fixed') on the Nevada tract.
Rebuilds the Nevada layers from raw data in area-relative mode, records the cut-offs each builder used (cap distances in metres, density and
fabric ranges), and writes reference/nevada_scale.json. Because Nevada's own layers were scaled with exactly these cut-offs, rebuilding Nevada in
fixed mode reproduces them (tests/test_fixed_scale_nevada.py). Needs the Terra-X Nevada project files (set FA)."""
import os, sys, json, warnings; warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import build_layers
FA = os.environ.get("FA", "."); OM = f"{FA}/OpenMineral"; PR = f"{FA}/Nevada Datasets/Processed Results"
OUT = os.environ.get("CAL_OUT", "/tmp/cal_nevada")
CAP = ["eng_volcanic_intrusive_proximity", "src_eocene_intrusive_proximity", "src_jurassic_cretaceous_intrusive_proximity", "arch_lithologic_contact_proximity",
       "src_basement_proximity", "trap_reactive_carbonate", "trap_impermeable_cap", "arch_fault_intersection_density"]
RANGE = ["arch_fault_density", "arch_structural_complexity", "arch_fluid_pathway_connectivity", "eng_eocene_belt_density"]

def nevada_cfg(scaling):
    cfg = json.load(open(os.path.join(HERE, "..", "references", "example_config.json"))); cfg["scaling"] = scaling; cfg["data_dir"] = OM
    cfg["tract"] = {"file": f"{PR}/Permissive Tract/Permissive_Tract.gpkg", "bounds": [273000, 4655000 - 725000, 273000 + 494000, 4655000]}; cfg["tract_is_permissive"] = True
    for k in ("arch_fault_density", "arch_fault_intersection_density", "arch_fluid_pathway_connectivity", "arch_structural_complexity"): cfg["layers"][k]["file"] = f"{PR}/Faults_NAD83.gpkg"
    cfg["layers"]["arch_lithologic_contact_proximity"] = {"file": "Geology_NAD83.shp", "field": "L_NAME", "regex": "."}
    for lid, c in cfg["layers"].items():
        if c.get("file") == "geology.shp": c["file"] = "Geology_NAD83.shp"
        if c.get("file") == "geochemistry.gpkg": c["file"] = "Geochemical_Data.gpkg"
        if "lithology" in c: c["lithology"]["file"] = "Geology_NAD83.shp"
        if "cover" in c: c["cover"]["file"] = "Geology_NAD83.shp"
        if c.get("k"): c["k"] = "K_26911.tif"
        if c.get("depth"): c["depth"] = "GB_depth_to_basement_surface.tif"
        if c.get("gravity"): c["gravity"] = "Gravity_Data.tif"
        if c.get("magnetics"): c["magnetics"] = "Magnetic Data.tif"
        if "intrusive" in c: c["intrusive"]["file"] = "Geology_NAD83.shp"
    return cfg

if __name__ == "__main__":
    only = sys.argv[1].split(",") if len(sys.argv) > 1 else CAP + RANGE
    cfg = nevada_cfg("relative"); cfg["layers"] = {k: v for k, v in cfg["layers"].items() if k in only}
    log = build_layers.build(cfg, os.path.join(OUT, "relative"))
    path = os.path.join(HERE, "..", "reference", "nevada_scale.json"); sc = json.load(open(path)) if os.path.exists(path) else {}
    for lid in only:
        r = log[lid]
        if r["result"] != "built": print("FAILED", lid, r["result"]); continue
        u = r["cutoffs_used"]
        if lid in CAP: sc[lid] = {"cap_m": u["cap_m"], "note": "distance cap in metres (Nevada 98th percentile, or farthest cell for intersections)"}
        else: sc[lid] = {"fixed": [u["lo"], u["hi"]], "note": "value range mapped to 0-1 (Nevada tract)"}
        print(lid, sc[lid])
    json.dump(sc, open(path, "w"), indent=1)
