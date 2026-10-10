"""Tests the newer builders on the raw OpenMineral data against the Nevada layers (read-only use of OpenMineral).
Set FA to the project folder."""
import os, sys, glob, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, geopandas as gpd, rasterio, rasterio.warp, rasterio.features
from scipy.stats import spearmanr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts")); import builders as B
FA = os.environ.get("FA", "."); OM = f"{FA}/OpenMineral"
LYR = f"{FA}/Analytical Results/Old Results/5b_Team_Handoff/02_Evidence_Layers"
def load(s):
    with rasterio.open(f"{LYR}/{s}.tif") as ds: a = ds.read(1).astype(float); a[a == ds.nodata] = np.nan; return a, ds
stems = [os.path.basename(p)[:-4] for p in sorted(glob.glob(f"{LYR}/*.tif"))]
tract = np.all(np.isfinite(np.stack([load(s)[0] for s in stems])), 0); _, ds = load(stems[0]); T = ds.transform; crs = ds.crs
grid = B.Grid((T.c, T.f + ds.height * T.e, T.c + ds.width * T.a, T.f))
res = []
def cmp(n, new, old):
    m = tract & np.isfinite(new) & np.isfinite(old); r = spearmanr(new[m], old[m]).correlation; res.append((n, round(r, 3))); print(f"{n:62s} Spearman {r:.3f}", flush=True)
def warp(path):
    with rasterio.open(path) as src:
        a = np.full((grid.H, grid.W), np.nan); rasterio.warp.reproject(rasterio.band(src, 1), a, src_transform=src.transform, src_crs=src.crs, dst_transform=grid.transform, dst_crs=crs, resampling=rasterio.warp.Resampling.bilinear, dst_nodata=np.nan, src_nodata=src.nodata)
    return np.where(a < -9000, np.nan, a)
# the OpenMineral fault file is corrupted: check the guard
try: B.clean_geoms(gpd.read_file(f"{OM}/Faults_NAD83.gpkg")); print("faults: no error (unexpected)")
except ValueError as e: print("faults guard OK:", str(e)[:80])
geo = gpd.read_file(f"{OM}/Geology_NAD83.shp"); geo["FMATN"] = geo["FMATN"].fillna("").astype(str)
codes = {c: i + 1 for i, c in enumerate(sorted(geo.FMATN.unique()))}
U = rasterio.features.rasterize([(g, codes[c]) for g, c in zip(geo.geometry, geo.FMATN) if g is not None and not g.is_empty], out_shape=(grid.H, grid.W), transform=grid.transform, fill=0, dtype="int32")
cmp("Potassium depletion (per-unit baseline)", B.potassium_depletion(grid, tract, warp(f"{OM}/K_26911.tif"), U), load("trap_argillic_k_depletion")[0])
cmp("Potassium depletion (ln K ~ ln Th, no lithology)", B.potassium_depletion_th(grid, tract, warp(f"{OM}/K_26911.tif"), warp(f"{OM}/Th_26911.tif")), load("trap_argillic_k_depletion")[0])
cover = B.polygons_to_mask(grid, list(geo[geo.FMATN.str.match(r"^(Q|QT)")].geometry))
cmp("Preservation (depth classes, Quaternary cover)", B.preservation_from_depth(grid, tract, warp(f"{OM}/GB_depth_to_basement_surface.tif"), cover), load("trap_preservation")[0])
gc = gpd.read_file(f"{OM}/Geochemical_Data.gpkg")
cmp("Magmatic fertility (Sr/Y, La/Y) vs Nevada layer", B.magmatic_fertility(grid, tract, gc), load("src_low_magmatic_fertility")[0])
cmp("Pathfinder ore-fluid footprint vs Nevada fluid chemistry", B.geochem_pathfinder_score(grid, tract, gc), load("trap_fluid_chemistry")[0])
ig = geo[geo.L_NAME.fillna("").str.contains("intrus|pluton|stock|dike|granit|rhyolit", case=False, regex=True)]
im = B.polygons_to_mask(grid, list(ig.geometry))
cmp("Buried plutons (gravity + magnetics)", B.buried_pluton_proximity(grid, tract, warp(f"{OM}/Gravity_Data.tif"), warp(f"{OM}/Magnetic Data.tif"), im & tract), load("eng_buried_pluton_proximity")[0])
pts = B.load_geochron(f"{OM}/ar40ar39_USGS_data_release_v4.csv", f"{OM}/kar_USGS_data_release_v4.csv", f"{OM}/samples_USGS_data_release_v4.csv", crs=str(crs))
w = pts[(pts.age >= 36) & (pts.age <= 42)]; print("USGS geochron points in 36-42 Ma:", len(w), "of", len(pts), "(Nevada report: 202 unique)")
nv = w.cx[grid.xmin:grid.xmax, grid.ymin:grid.ymax]; print("inside the Nevada grid:", len(nv))
cmp("Belt density from raw USGS tables", B.age_window_density(grid, tract, nv.geometry.x.values, nv.geometry.y.values, nv.age.values, 36, 42), load("eng_eocene_belt_density")[0])
pd.DataFrame(res, columns=["layer", "spearman_vs_original"]).to_csv(os.path.join(os.path.dirname(__file__), "builders_openmineral_results.csv"), index=False)
