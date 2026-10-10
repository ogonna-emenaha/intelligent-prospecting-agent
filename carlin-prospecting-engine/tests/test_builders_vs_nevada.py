"""Rebuild Nevada layers with the generic builders and compare with the original layers (Spearman on tract cells).
Needs the Terra-X Nevada project files. Set FA to the project folder."""
import os, sys, glob, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from scipy.stats import spearmanr
from scipy import ndimage as ndi
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts")); import builders as B
FA = os.environ.get("FA", ".")
NV = f"{FA}/Nevada Datasets/Processed Results"; LYR = f"{FA}/Analytical Results/Old Results/5b_Team_Handoff/02_Evidence_Layers"
def load(s):
    with rasterio.open(f"{LYR}/{s}.tif") as ds:
        a = ds.read(1).astype(float); a[a == ds.nodata] = np.nan; return a, ds
stems = [os.path.basename(p)[:-4] for p in sorted(glob.glob(f"{LYR}/*.tif"))]
tract = np.all(np.isfinite(np.stack([load(s)[0] for s in stems])), 0); _, ds = load(stems[0]); T = ds.transform
grid = B.Grid((T.c, T.f + ds.height * T.e, T.c + ds.width * T.a, T.f)); assert (grid.H, grid.W) == tract.shape
res = []
def cmp(name, new, old):
    m = tract & np.isfinite(new) & np.isfinite(old); r = spearmanr(new[m], old[m]).correlation; res.append((name, round(r, 3))); print(f"{name:58s} Spearman {r:.3f}", flush=True)
geo = gpd.read_file(f"{NV}/Geology_Fixed.gpkg")
ig = geo[geo.L_NAME.fillna("").str.contains("intrus|pluton|stock|dike|granit|rhyolit", case=False, regex=True)]
cmp("Proximity to volcanic and intrusive rocks", B.proximity_score(B.polygons_to_mask(grid, ig.geometry), tract), load("eng_volcanic_intrusive_proximity")[0])
pts = gpd.read_file(f"{NV}/Energy Analysis/eocene_dated_points_26911.gpkg"); age = pts["age"].astype(float).values
print("dated points", len(pts), "age range", np.nanmin(age), np.nanmax(age), flush=True)
x = pts.geometry.x.values; y = pts.geometry.y.values
cmp("Intrusive belt density, 36-42 Ma window", B.age_window_density(grid, tract, x, y, age, 36, 42), load("eng_eocene_belt_density")[0])
cmp("Proximity to intrusions, 36-42 Ma window", B.age_window_proximity(grid, tract, x, y, age, 36, 42), load("src_eocene_intrusive_proximity")[0])
fl = gpd.read_file(f"{NV}/Faults_NAD83.gpkg"); print(fl.F_TYPE.value_counts().to_dict(), flush=True)
def wt(t):
    t = str(t).lower()
    if "thrust" in t: return 0.5 if "known" in t else 0.35
    if "inferred" in t: return 0.7
    if "conceal" in t: return 0.4
    return 1.0
w = fl.F_TYPE.map(wt).values
cmp("Fault density (weighted line density, 10 km)", B.fault_density(grid, tract, list(fl.geometry), w), load("arch_fault_density")[0])
ix = gpd.read_file(f"{NV}/Pathway Analysis/Layer2_FaultIntersectionDensity/fault_intersections.gpkg")
m = B.points_to_mask(grid, ix.geometry.x.values, ix.geometry.y.values)
d = ndi.distance_transform_edt(~m) * grid.cell; cmp("Fault intersection proximity (supplied intersections)", np.where(tract, 1 - d / np.nanmax(d[tract]), np.nan), load("arch_fault_intersection_density")[0])
sc = B.structural_complexity(grid, tract, list(fl.geometry))
if sc is not None: cmp("Structural complexity (1 - R, 10 km)", sc, load("arch_structural_complexity")[0])
pd.DataFrame(res, columns=["layer", "spearman_vs_original"]).to_csv(os.path.join(os.path.dirname(__file__), "builder_vs_nevada_results.csv"), index=False)

# ---- second batch: contact, basement, older intrusions, gravity layers ----
import re
geo["FMATN"] = geo["FMATN"].fillna("").astype(str)
bs = geo[geo.FMATN.str.startswith("pC") & ~geo.FMATN.isin(["pCg"])]
cmp("Distance from basement outcrops (pC units)", B.basement_distance(grid, tract, bs.geometry), load("src_basement_proximity")[0])
jk = geo[geo.FMATN.str.match(r"^[JK]") & geo.L_NAME.fillna("").str.contains("granit|pluton|diorit|intrusiv|stock|quartz monzonit", case=False, regex=True)]
cmp("Proximity to older (Jurassic-Cretaceous) intrusions", B.proximity_score(B.polygons_to_mask(grid, jk.geometry), tract), load("src_jurassic_cretaceous_intrusive_proximity")[0])
cmp("Proximity to lithologic contacts (unit boundaries)", B.proximity_score(B.boundary_mask(grid, [g for g in geo.geometry if g is not None]), tract), load("arch_lithologic_contact_proximity")[0])
import rasterio.warp
GR = f"{FA}/Nevada Datasets/GB_geophysics_grids/GB_geophysics_grids/GB_iso_grav_anom.tif"
with rasterio.open(GR) as src:
    g = np.full(tract.shape, np.nan); rasterio.warp.reproject(rasterio.band(src, 1), g, src_transform=src.transform, src_crs=src.crs, dst_transform=T, dst_crs=ds.crs, resampling=rasterio.warp.Resampling.bilinear, dst_nodata=np.nan)
g = np.where(g < -9000, np.nan, g)
cmp("Mantle upwelling (gravity upward-continued 8 km)", B.upward_continuation(g, tract, 8000), load("eng_mantle_upwelling")[0])
cmp("Crustal gradient (log total horizontal gradient)", B.total_horizontal_gradient(g, tract), load("src_crustal_gradient_fav")[0])
pd.DataFrame(res, columns=["layer", "spearman_vs_original"]).to_csv(os.path.join(os.path.dirname(__file__), "builder_vs_nevada_results.csv"), index=False)
