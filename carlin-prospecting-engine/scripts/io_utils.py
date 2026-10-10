import os, glob, json, numpy as np, rasterio
def read_layers(folder, wanted_ids):
    """Read GeoTIFF layers named <layer_id>.tif from folder. All must share one grid. Returns dict, profile."""
    L = {}; prof = None
    for s in wanted_ids:
        p = os.path.join(folder, s + ".tif")
        if not os.path.exists(p): continue
        with rasterio.open(p) as ds:
            a = ds.read(1).astype(float)
            if ds.nodata is not None: a[a == ds.nodata] = np.nan
            a[a < -9000] = np.nan
            if prof is None: prof = dict(transform=ds.transform, crs=ds.crs, shape=a.shape)
            elif a.shape != prof["shape"] or ds.transform != prof["transform"]: raise ValueError(f"{s}.tif is not on the same grid as the other layers")
            L[s] = a
    return L, prof
def write_raster(path, arr, prof, nodata=-9999.0, dtype="float32"):
    out = np.where(np.isfinite(arr), arr, nodata).astype(dtype)
    with rasterio.open(path, "w", driver="GTiff", height=out.shape[0], width=out.shape[1], count=1, dtype=dtype, crs=prof["crs"], transform=prof["transform"], nodata=nodata) as ds: ds.write(out, 1)

def deposits_xy(path, layers_dir):
    """Read deposit points (gpkg/shp/geojson/csv) and return x,y in the CRS of the layers in layers_dir."""
    import geopandas as gpd, pandas as pd, glob
    with rasterio.open(sorted(glob.glob(os.path.join(layers_dir, "*.tif")))[0]) as ds: crs = ds.crs
    if path.lower().endswith(".csv"):
        d = pd.read_csv(path); cols = {c.lower(): c for c in d.columns}
        xc = next(cols[k] for k in ("x", "easting", "lon", "longitude") if k in cols); yc = next(cols[k] for k in ("y", "northing", "lat", "latitude") if k in cols)
        geo = gpd.GeoSeries(gpd.points_from_xy(d[xc], d[yc]), crs="EPSG:4326" if xc.lower() in ("lon", "longitude") else crs)
    else:
        geo = gpd.read_file(path).geometry
    geo = geo.to_crs(crs); return np.c_[geo.x.values, geo.y.values]
