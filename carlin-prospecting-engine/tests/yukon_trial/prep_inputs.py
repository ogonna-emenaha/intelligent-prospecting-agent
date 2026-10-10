"""Prepare Yukon inputs for the new engine WITHOUT touching the OpenMineral folder (read-only). Masks of the geology-derived layers are recovered from
the earlier run's layers (cells scoring exactly 1 for proximity layers, 0 for the basement-distance layer) because the unit-selection files were not saved."""
import geopandas as gpd, rasterio, numpy as np, rasterio.features, shapely, os
from shapely.geometry import shape
Y=os.environ.get("YUKON","Yukon_Prospectivity"); R=f"{Y}/Carlin_engine_run_2026-10-08"; D=f"{R}/inputs_derived"; OUT="inputs"; os.makedirs(OUT,exist_ok=True)
def mask_poly(layer, inside_value):
    with rasterio.open(f"{R}/layers/{layer}.tif") as ds:
        a=ds.read(1); m=(a>-9000)&(np.abs(a-inside_value)<1e-6); T=ds.transform
    gs=[shape(g) for g,v in rasterio.features.shapes(m.astype("uint8"),mask=m,transform=T)]
    gpd.GeoDataFrame({"SEL":["x"]*len(gs)},geometry=gs,crs=3578).to_file(f"{OUT}/{layer}.gpkg",driver="GPKG"); print(layer,len(gs),int(m.sum()),"cells")
for l,v in [("trap_reactive_carbonate",1),("trap_impermeable_cap",1),("eng_volcanic_intrusive_proximity",1),("src_jurassic_cretaceous_intrusive_proximity",1),("src_basement_proximity",0)]: mask_poly(l,v)
f=gpd.read_file(f"{D}/faults_classed.gpkg"); raw=gpd.read_file(f"{Y}/Faults.shp"); assert len(f)==len(raw); f["SCALE"]=raw["SCALE"].values; f.to_file(f"{OUT}/faults_classed_scale.gpkg",driver="GPKG")
c=gpd.read_file(f"{D}/reactive_cap_contacts.gpkg"); c["SEL"]="x"; c.to_file(f"{OUT}/reactive_cap_contacts.gpkg",driver="GPKG")
w=gpd.read_file(f"{D}/window_igneous_points.gpkg"); w.to_file(f"{OUT}/window_igneous_points.gpkg",driver="GPKG")
print("scale classes",f.SCALE.value_counts().to_dict())
