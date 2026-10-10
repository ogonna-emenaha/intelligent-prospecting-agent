"""Smoke test of the whole chain on Nevada raw data (OpenMineral data read-only; faults from the project's working copy because the
OpenMineral fault file is corrupted). build layers -> zero-shot -> compare with the stored-layer reference.
NOT an independent validation: the Nevada deposits are in the reference."""
import os, sys, json, warnings; warnings.filterwarnings("ignore")
import numpy as np, rasterio
from scipy.stats import spearmanr
HERE = os.path.dirname(__file__); sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import build_layers, run_zero_shot
FA = os.environ.get("FA", "."); OM = f"{FA}/OpenMineral"; PR = f"{FA}/Nevada Datasets/Processed Results"
OUT = os.environ.get("E2E_OUT", "/tmp/e2e_nevada"); os.makedirs(OUT, exist_ok=True)
cfg = json.load(open(os.path.join(HERE, "..", "references", "example_config.json")))
cfg["data_dir"] = OM; cfg["tract_is_permissive"] = True; cfg["tract_basis"] = "Nevada Carlin permissive tract (Terra-X)"
cfg["tract"] = {"file": f"{PR}/Permissive Tract/Permissive_Tract.gpkg", "bounds": [273000, 4655000 - 725000, 273000 + 494000, 4655000]}
for k in ("arch_fault_density", "arch_fault_intersection_density", "arch_fluid_pathway_connectivity", "arch_structural_complexity"): cfg["layers"][k]["file"] = f"{PR}/Faults_NAD83.gpkg"
cfg["layers"]["arch_lithologic_contact_proximity"] = {"file": "Geology_NAD83.shp", "field": "L_NAME", "regex": "."}
for lid, c in cfg["layers"].items():
    for key in ("file",):
        if c.get(key) == "geology.shp": c[key] = "Geology_NAD83.shp"
        if c.get(key) == "geochemistry.gpkg": c[key] = "Geochemical_Data.gpkg"
    if "lithology" in c: c["lithology"]["file"] = "Geology_NAD83.shp"
    if "cover" in c: c["cover"]["file"] = "Geology_NAD83.shp"
    if c.get("k"): c["k"] = "K_26911.tif"
    if c.get("depth"): c["depth"] = "GB_depth_to_basement_surface.tif"
    if c.get("gravity"): c["gravity"] = "Gravity_Data.tif"
    if c.get("magnetics"): c["magnetics"] = "Magnetic Data.tif"
    if "intrusive" in c: c["intrusive"]["file"] = "Geology_NAD83.shp"
log = build_layers.build(cfg, f"{OUT}/layers"); print(json.dumps({k: v.get("result") + " | " + str(v.get("proxy_status")) for k, v in log.items() if k != "_meta"}, indent=1))
rep = run_zero_shot.run(f"{OUT}/layers", f"{OUT}/zero_shot", n_bags=10)
for k in ("verdict", "n_cells", "share_outside_reference_domain", "parts_missing", "layers_missing"): print(k, ":", rep.get(k))
print("proxy layers:", [(p["layer"][:40], p["status"]) for p in rep["proxy_layers"]]); print("layers to check:", [(x["layer"], x["district_median"], x["nevada_p5"], x["nevada_p95"]) for x in rep["layers_to_check"]])
LYR = f"{FA}/Analytical Results/Old Results/5b_Team_Handoff/02_Evidence_Layers"
for lid, st in log.items():
    if lid == "_meta": continue
    if st["result"] != "built": continue
    with rasterio.open(f"{OUT}/layers/{lid}.tif") as a, rasterio.open(f"{LYR}/{lid}.tif") as b: x = a.read(1); y = b.read(1)
    m = (x > -9000) & (y > -9000) & np.isfinite(x) & np.isfinite(y); print(f"{lid:46s} rank r={spearmanr(x[m], y[m]).correlation:6.3f}")
REF = os.path.join(HERE, "..", "reference")
with rasterio.open(f"{OUT}/zero_shot/carlin_similarity_percentile.tif") as a, rasterio.open(f"{REF}/reference_magmatic_prospectivity_engine.tif") as b: x = a.read(1); y = b.read(1)
m = (x > -9000) & (y > -9000); print("Spearman, rebuilt-from-raw zero-shot vs stored-layer reference:", round(spearmanr(x[m], y[m]).correlation, 3))
