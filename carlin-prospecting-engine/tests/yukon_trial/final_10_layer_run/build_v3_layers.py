"""Version 3 layers (whole Yukon, 1 km, EPSG:3578). Builds the two intrusion layers with Nevada fixed scaling from 88-102 Ma intrusive polygons,
rebuilds the As-Sb-Au ore-fluid layer over the whole Yukon, copies the 7 structural and host-rock layers from the earlier build, runs nothing else."""
import sys, shutil, numpy as np, pandas as pd, geopandas as gpd, rasterio
from scipy import ndimage as ndi
sys.path.insert(0, "ENGINE_SCRIPTS_DIR"); import builders as B          # Carlin engine scripts folder
from io_utils import write_raster
with rasterio.open("out/layers8/trap_reactive_carbonate.tif") as d: ref = d.read(1); T = d.transform; prof = dict(transform=T, crs=d.crs, shape=ref.shape)
tract = ref > -9000; H, W = ref.shape; grid = B.Grid((T.c, T.f - H * 1000, T.c + W * 1000, T.f), 1000, "EPSG:3578")
g = gpd.read_file("Bedrock_Geology.shp"); mid = (g.AGE_MAX_MA + g.AGE_MIN_MA) / 2
s = g[g.ROCK_CLASS.str.contains("plut|intrus", case=False, na=False) & mid.between(88, 102)].copy(); s["MID_AGE"] = mid[s.index]   # 688 polygons
s[["GP_SUITE", "UNIT_1M", "AGE_MAX_MA", "AGE_MIN_MA", "MID_AGE", "geometry"]].to_file("v3/source_intrusions_88_102Ma.gpkg")
# Source: distance to nearest selected polygon, capped at 275 km, score 1 - d/275 km
d = ndi.distance_transform_edt(~B.polygons_to_mask(grid, list(s.geometry))) * 1000.0
write_raster("v3/layers/src_jurassic_cretaceous_intrusive_proximity.tif", np.where(tract, 1 - np.minimum(d, 275000.0) / 275000.0, np.nan), prof)
# Energy: one representative point per polygon, 140 km quartic kernel, samples per km2, fixed range 0 to 0.00352956 (clipped above)
rp = s.representative_point()
bd = B.age_window_density(grid, tract, rp.x.values, rp.y.values, s.MID_AGE.values, 88, 102, 140000, (0.0, 0.00352956))
write_raster("v3/layers/eng_eocene_belt_density.tif", bd, prof)
for l in ["arch_fault_density", "arch_fault_intersection_density", "arch_fluid_pathway_connectivity", "arch_structural_complexity",
          "arch_lithologic_contact_proximity", "trap_reactive_carbonate", "trap_impermeable_cap"]: shutil.copy(f"out/layers8/{l}.tif", "v3/layers/")
# Trap: ore-fluid footprint, log(geomean As, Sb, Au / geomean Th, Hf, Sc), IDW k=12 over the whole Yukon
rg = gpd.read_file("prep/rgs_clean.gpkg")
write_raster("v3/layers/trap_fluid_chemistry.tif", B.geochem_pathfinder_score(grid, tract, rg, ["AS_INA_PPM", "SB_INA_PPM", "AU_INA_PPB"], ["TH_INA_PPM", "HF_INA_PPM", "SC_INA_PPM"]), prof)
# then:  python run_zero_shot.py --layers v3/layers --out v3/zero_shot_magmatic_10 --deposits prep/validation_carlin_ygs.gpkg
