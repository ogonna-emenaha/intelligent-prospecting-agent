"""Generic, region-agnostic evidence-layer builders. Each returns a float array on a common grid (0-1, nan outside tract).
Conventions copied from the Nevada build (1 km cells, scores rescaled 0-1, relative to the study tract)."""
import numpy as np
from scipy import ndimage as ndi
from scipy.signal import fftconvolve

def _conv(a, k):
    return np.maximum(fftconvolve(a, k, mode='same'), 0.0)


class Grid:
    """Common analysis grid. bounds=(xmin,ymin,xmax,ymax) in a projected CRS in metres."""
    def __init__(self, bounds, cell=1000.0, crs="EPSG:26911"):
        self.xmin, self.ymin, self.xmax, self.ymax = bounds; self.cell = cell; self.crs = crs
        self.W = int(np.ceil((self.xmax - self.xmin) / cell)); self.H = int(np.ceil((self.ymax - self.ymin) / cell))
    @property
    def transform(self):
        import rasterio
        return rasterio.transform.from_origin(self.xmin, self.ymin + self.H * self.cell, self.cell, self.cell)
    def rc(self, x, y):
        c = np.floor((np.asarray(x) - self.xmin) / self.cell).astype(int)
        r = np.floor((self.ymin + self.H * self.cell - np.asarray(y)) / self.cell).astype(int)
        return r, c

# ---- scaling ----------------------------------------------------------------------------------------------------------
# Two modes. Area-relative (legacy): each layer is rescaled 0-1 using the min/max or percentile cut-off of the area it is run on.
# Fixed (default in build_layers): the same physical cut-offs are used in every region (calibrated once on the Nevada tract and stored in
# reference/nevada_scale.json), so a score of 0.8 means the same distance, density or fabric disorder in every district.
# LAST records the cut-offs used by the most recent builder call (read by build_layers and by the Nevada calibration).
LAST = {}

def minmax(a, tract, fixed=None):
    v = a[tract]; lo, hi = (np.nanmin(v), np.nanmax(v)) if fixed is None else (float(fixed[0]), float(fixed[1]))
    LAST.update(lo=float(lo), hi=float(hi)); out = (a - lo) / (hi - lo + 1e-12)
    if fixed is not None: out = np.clip(out, 0.0, 1.0)
    return np.where(tract, out, np.nan)

def winsor_minmax(a, tract, lo_q=2, hi_q=98, fixed=None):
    if fixed is not None: return minmax(a, tract, fixed)
    v = a[tract]; lo, hi = np.nanpercentile(v, [lo_q, hi_q]); a = np.clip(a, lo, hi); return minmax(a, tract)

def proximity_score(mask, tract, winsor_q=98, cell=1000.0, cap_m=None):
    """Euclidean distance to the nearest True cell, capped at cap_m, then inverted: score = 1 - min(d, cap)/cap.
    cap_m=None (area-relative): the cap is the winsor_q-th percentile of the distances inside the tract. cap_m given (fixed scale): the
    Nevada-calibrated cap in metres. (Nevada convention for volcanic/intrusive, intrusion, contact, basement layers.)"""
    d = ndi.distance_transform_edt(~mask) * cell
    p = np.nanpercentile(d[tract], winsor_q) if cap_m is None else float(cap_m)
    LAST.update(cap_m=float(p))
    s = 1 - np.minimum(d, p) / p
    return np.where(tract, s, np.nan)

def points_to_mask(grid, x, y):
    m = np.zeros((grid.H, grid.W), bool); r, c = grid.rc(x, y)
    ok = (r >= 0) & (r < grid.H) & (c >= 0) & (c < grid.W); m[r[ok], c[ok]] = True; return m

def lines_to_mask(grid, geoms, step=None):
    """Rasterise line geometries by densifying vertices every half cell."""
    step = step or grid.cell / 2; m = np.zeros((grid.H, grid.W), bool)
    for g in geoms:
        if g is None or g.is_empty: continue
        parts = getattr(g, "geoms", [g])
        for p in parts:
            n = max(int(p.length // step), 1) + 1
            pts = [p.interpolate(i / (n - 1), normalized=True) for i in range(n)] if n > 1 else [p.interpolate(0)]
            r, c = grid.rc([q.x for q in pts], [q.y for q in pts]); ok = (r >= 0) & (r < grid.H) & (c >= 0) & (c < grid.W); m[r[ok], c[ok]] = True
    return m

def polygons_to_mask(grid, geoms):
    import rasterio.features
    shapes = [(g, 1) for g in geoms if g is not None and not g.is_empty]
    if not shapes: return np.zeros((grid.H, grid.W), bool)
    return rasterio.features.rasterize(shapes, out_shape=(grid.H, grid.W), transform=grid.transform, fill=0, dtype="uint8").astype(bool)

def boundary_mask(grid, polys):
    return lines_to_mask(grid, [p.boundary for p in polys if p is not None and not p.is_empty])

def quartic_kde(points_mask_counts, radius_m, cell=1000.0):
    """Quartic (biweight) kernel density of a count grid with search radius radius_m. Returned in samples per km2 (kernel integrates to 1),
    so the value does not depend on the cell size."""
    R = int(round(radius_m / cell)); y, x = np.mgrid[-R:R + 1, -R:R + 1]; d = np.hypot(x, y) / R
    k = np.where(d <= 1, (1 - d ** 2) ** 2, 0.0); k = k / k.sum() / (cell / 1000.0) ** 2; return _conv(points_mask_counts.astype(float), k)

def count_grid(grid, x, y):
    c = np.zeros((grid.H, grid.W)); r, cc = grid.rc(x, y); ok = (r >= 0) & (r < grid.H) & (cc >= 0) & (cc < grid.W)
    np.add.at(c, (r[ok], cc[ok]), 1); return c

# ---- named layer builders ----
def age_window_density(grid, tract, x, y, age, age_min, age_max, radius_m=140000, fixed=None):
    """Intrusive belt density: quartic KDE of dated igneous samples inside [age_min, age_max] Ma (user-set window).
    NOTE: the Nevada report documents a 15 km radius, but the delivered Nevada raster is reproduced (Spearman 0.995) only by a
    radius near 140 km, so 140 km is the default. The layer therefore measures regional magmatic-belt presence, not local density."""
    sel = (age >= age_min) & (age <= age_max)
    return minmax(quartic_kde(count_grid(grid, x[sel], y[sel]), radius_m, grid.cell), tract, fixed)

def age_window_proximity(grid, tract, x, y, age, age_min, age_max, cap_m=None):
    sel = (age >= age_min) & (age <= age_max)
    return proximity_score(points_to_mask(grid, x[sel], y[sel]), tract, cap_m=cap_m)

def fault_density(grid, tract, geoms, weights=None, radius_m=10000, fixed=None):
    """Weighted fault length per unit area in a circular neighbourhood (approximately km of fault per km2). weights: optional per-geometry
    weight (Nevada: known 1.0, inferred 0.7, ...)."""
    w = np.ones(len(geoms)) if weights is None else np.asarray(weights, float); dens = np.zeros((grid.H, grid.W))
    R = int(round(radius_m / grid.cell)); y, x = np.mgrid[-R:R + 1, -R:R + 1]; k = (np.hypot(x, y) <= R).astype(float)
    for wi in np.unique(w):
        m = lines_to_mask(grid, [g for g, ww in zip(geoms, w) if ww == wi]); dens += wi * _conv(m.astype(float), k)
    dens = dens / (k.sum() * grid.cell / 1000.0)
    return minmax(dens, tract, fixed)

def generalise_lines(geoms, scales, target_scale=250, tol_m=125.0, min_len_m=2000.0):
    """Bring faults mapped at different detail to one common mapping scale. scales: per-feature map-scale denominator in thousands (for example
    250 for 1:250,000, 50 for 1:50,000). Features mapped in more detail than target_scale are simplified (tolerance tol_m, 0.5 mm at 1:250,000) and
    dropped if shorter than min_len_m, which is roughly what a 1:250,000 compilation would have kept. Returns (geoms, kept_mask)."""
    out, keep = [], []
    for g, s in zip(geoms, scales):
        try: finer = float(s) < target_scale
        except Exception: finer = False
        if finer:
            g = g.simplify(tol_m, preserve_topology=False)
            if g.is_empty or g.length < min_len_m: out.append(g); keep.append(False); continue
        out.append(g); keep.append(True)
    return out, np.array(keep, bool)

def fault_intersections(geoms, geoms_b=None):
    """Intersection points between fault geometries. If geoms_b is given, only A x B intersections (Nevada: high-angle conduit
    faults x thrust faults, 5,181 points); otherwise all pairs (approximate)."""
    from shapely.geometry import Point
    from shapely.strtree import STRtree
    pts = []; other = geoms if geoms_b is None else geoms_b; tree = STRtree(other)
    for i, g in enumerate(geoms):
        for j in tree.query(g):
            if geoms_b is None and j <= i: continue
            inter = g.intersection(other[j])
            if inter.is_empty: continue
            for p in getattr(inter, "geoms", [inter]):
                if isinstance(p, Point): pts.append((p.x, p.y))
    return np.array(pts)

def fault_intersection_proximity(grid, tract, geoms, geoms_b=None, cap_m=None):
    """Nevada convention (validated AUC 0.786 vs 0.651 for KDE): distance to nearest intersection, inverted 1 - d/dmax.
    dmax is the largest distance inside the tract (area-relative) or the fixed Nevada cap cap_m."""
    pts = fault_intersections(geoms, geoms_b); m = points_to_mask(grid, pts[:, 0], pts[:, 1])
    d = ndi.distance_transform_edt(~m) * grid.cell; dm = np.nanmax(d[tract]) if cap_m is None else float(cap_m); LAST.update(cap_m=float(dm))
    return np.where(tract, np.clip(1 - d / dm, 0.0, 1.0), np.nan)

def structural_complexity(grid, tract, geoms, radius_m=10000, step=100.0, fixed=None):
    """1 - mean resultant length R of local fault strike (axial data, 10 km radius). High = multi-directional fabric."""
    C = np.zeros((grid.H, grid.W)); S = np.zeros((grid.H, grid.W)); N = np.zeros((grid.H, grid.W))
    for g in geoms:
        for p in getattr(g, "geoms", [g]):
            n = max(int(p.length // step), 1) + 1
            pts = [p.interpolate(i / (n - 1), normalized=True) for i in range(n)] if n > 1 else [p.interpolate(0)]
            xs = np.array([q.x for q in pts]); ys = np.array([q.y for q in pts])
            if len(xs) < 2: continue
            th = np.arctan2(np.gradient(ys), np.gradient(xs)); r, c = grid.rc(xs, ys); ok = (r >= 0) & (r < grid.H) & (c >= 0) & (c < grid.W)
            np.add.at(C, (r[ok], c[ok]), np.cos(2 * th[ok])); np.add.at(S, (r[ok], c[ok]), np.sin(2 * th[ok])); np.add.at(N, (r[ok], c[ok]), 1)
    R = int(round(radius_m / grid.cell)); y, x = np.mgrid[-R:R + 1, -R:R + 1]; k = (np.hypot(x, y) <= R).astype(float)
    Cs, Ss, Ns = (fftconvolve(a, k, mode='same') for a in (C, S, N))
    Rbar = np.where(Ns > 0, np.hypot(Cs, Ss) / np.maximum(Ns, 1), np.nan)
    return minmax(1 - Rbar, tract & np.isfinite(Rbar), fixed) if np.isfinite(Rbar[tract]).any() else None

def upward_continuation(grav, tract, height_m, cell=1000.0):
    """Upward continue a gravity grid by height_m (FFT), then rescale. Nevada 'mantle upwelling' used 8000 m."""
    g = np.where(np.isfinite(grav), grav, np.nanmean(grav)); ny, nx = g.shape
    kx = np.fft.fftfreq(nx, cell) * 2 * np.pi; ky = np.fft.fftfreq(ny, cell) * 2 * np.pi
    K = np.hypot(*np.meshgrid(kx, ky)); out = np.real(np.fft.ifft2(np.fft.fft2(g) * np.exp(-K * height_m)))
    return minmax(out, tract)

def total_horizontal_gradient(grav, tract, cell=1000.0):
    """Favourability = 1 - rescaled log THG. THG = sqrt(dg/dx^2 + dg/dy^2), log10, 2nd-98th percentile winsorised, rescaled 0-1 (Nevada 'crustal gradient')."""
    gy, gx = np.gradient(np.where(np.isfinite(grav), grav, np.nanmean(grav)), cell)
    thg = np.log10(np.hypot(gx, gy) + 1e-12); return 1 - winsor_minmax(thg, tract, 2, 98)   # favourability = flat (low-gradient) crust, as in Nevada

def basement_distance(grid, tract, basement_polys, cap_m=None):
    """Nevada layer 'basement proximity' is stored as distance FROM basement outcrop (higher = farther): deposits sit away from exposures."""
    return 1 - proximity_score(polygons_to_mask(grid, basement_polys), tract, cap_m=cap_m)

# ================= additional builders (tested against the Nevada layers where noted) =================
def clean_geoms(gdf, kind="lines"):
    """Drop None, empty and non-finite geometries. Raises a clear error if the file is corrupted (for example every coordinate infinite)."""
    import shapely
    g = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty]
    c = shapely.get_coordinates(g.geometry.values)
    if len(c) and not np.isfinite(c).any():
        raise ValueError("every coordinate in this file is infinite or NaN; the file is corrupted (often a failed reprojection or export). Re-export it from the source.")
    bb = g.geometry.bounds.to_numpy(); return g[np.isfinite(bb).all(1)]

def fault_network_connectivity(grid, tract, geoms, radius_km=5, fixed=None):
    """Approximate (0.56 vs the Nevada layer; the Nevada definition was not recovered): total length of the connected fault network
    each fault belongs to, averaged within radius_km, log-scaled."""
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    import rasterio.features
    from shapely.strtree import STRtree
    tree = STRtree(geoms); pr = tree.query(geoms, predicate="intersects"); m = pr[0] < pr[1]
    A = sp.coo_matrix((np.ones(m.sum()), (pr[0][m], pr[1][m])), shape=(len(geoms), len(geoms))); _, lab = connected_components(A, directed=False)
    L = np.array([g.length for g in geoms]); per = np.bincount(lab, weights=L)[lab]; order = np.argsort(per)
    res = rasterio.features.rasterize([(geoms[i], float(per[i])) for i in order], out_shape=(grid.H, grid.W), transform=grid.transform, fill=0, all_touched=True, dtype="float64")
    return minmax(np.log1p(ndi.uniform_filter(res, size=2 * int(radius_km) + 1)), tract, fixed)

def potassium_depletion(grid, tract, K, unit_raster=None):
    """Alteration proxy: potassium lower than expected for the rock type. With a lithology-unit raster (integer codes) the baseline is the
    unit's median ln K (0.66 vs the Nevada layer); without it, a regional ln K ~ ln Th regression residual is used (0.31)."""
    ok = tract & np.isfinite(K) & (K > 0); lk = np.where(ok, np.log(np.where(K > 0, K, 1)), np.nan)
    if unit_raster is not None:
        import pandas as pd
        ok &= unit_raster > 0; base = pd.Series(lk[ok]).groupby(unit_raster[ok]).transform("median").values; r = np.full(K.shape, np.nan); r[ok] = lk[ok] - base
    else:
        raise ValueError("pass a lithology unit raster, or use potassium_depletion_th")
    return winsor_minmax(np.where(ok, -r, np.nan), tract & ok, 1, 99)

def potassium_depletion_th(grid, tract, K, Th):
    ok = tract & np.isfinite(K) & np.isfinite(Th) & (K > 0) & (Th > 0); x = np.log(Th[ok]); y = np.log(K[ok]); b, a0 = np.polyfit(x, y, 1)
    r = np.full(K.shape, np.nan); r[ok] = y - (a0 + b * x); return winsor_minmax(np.where(ok, -r, np.nan), tract & ok, 1, 99)

def _num(s):
    import pandas as pd
    s = s.astype(str).str.strip(); cens = s.str.startswith("<"); v = pd.to_numeric(s.str.replace("<", "").str.replace(">", ""), errors="coerce"); return v.where(~cens, v / 2)

def idw_grid(grid, tract, x, y, v, k=12, power=2):
    from scipy.spatial import cKDTree
    t = cKDTree(np.c_[x, y]); gx = grid.xmin + (np.arange(grid.W) + .5) * grid.cell; gy = grid.ymin + grid.H * grid.cell - (np.arange(grid.H) + .5) * grid.cell
    GX, GY = np.meshgrid(gx, gy); q = np.c_[GX[tract], GY[tract]]; d, i = t.query(q, k=k); w = 1 / np.maximum(d, 1.0) ** power
    out = np.full((grid.H, grid.W), np.nan); out[tract] = (w * v[i]).sum(1) / w.sum(1); return out

def geochem_pathfinder_score(grid, tract, gdf, pathfinders=("As_ppm", "Sb_ppm", "Tl_ppm"), immobiles=("Al_pct", "Ti_pct", "Th_ppm")):
    """Ore-fluid footprint from pathfinder geochemistry: log of geometric-mean pathfinders over geometric-mean immobile elements, winsorised 1-99,
    IDW (k=12, power 2) to the grid. Element names are the user's column names. Different from Nevada's Au/Ag proxy (0.29 agreement), so treat as a proxy substitution."""
    cols = list(pathfinders) + list(immobiles); d = gdf.copy()
    for c in cols: d[c] = _num(d[c])
    m = np.all([d[c] > 0 for c in cols], axis=0); gm = lambda cs: np.exp(np.mean([np.log(d[c][m]) for c in cs], axis=0))
    z = np.log(gm(pathfinders) / gm(immobiles)); lo, hi = np.percentile(z, [1, 99]); z = np.clip(z, lo, hi)
    return minmax(idw_grid(grid, tract, d.geometry.x[m].values, d.geometry.y[m].values, np.asarray(z)), tract)

def gold_silver_footprint(grid, tract, gdf, au="Au_ppm", ag="Ag_ppm", base=("Cu_ppm", "Pb_ppm", "Zn_ppm"), ratio=10.0):
    """Nevada proxy: high Au/Ag (above `ratio`) together with absence of base metals. Needs gold assays."""
    d = gdf.copy()
    for c in [au, ag, *base]: d[c] = _num(d[c])
    m = (d[au] > 0) & (d[ag] > 0); r = np.log10(d[au][m] / d[ag][m]) - np.log10(ratio); bm = np.mean([np.log(d[c][m].clip(lower=1e-3)) for c in base], axis=0)
    s = 1 / (1 + np.exp(-r)) * (1 - minmax_arr(bm)); return minmax(idw_grid(grid, tract, d.geometry.x[m].values, d.geometry.y[m].values, np.asarray(s)), tract)

def minmax_arr(v):
    v = np.asarray(v, float); return (v - v.min()) / (v.max() - v.min() + 1e-12)

def magmatic_fertility(grid, tract, gdf, sr="Sr_ppm", y="Y_ppm", la="La_ppm"):
    """Magmatic fertility index 0.5 ln(Sr/Y) + 0.5 ln(La/Y), winsorised 1-99, IDW to the grid. High = more fertile (adakitic).
    In Nevada, deposits sit where this index is LOW (AUC 0.29), which the model learns; the layer itself is not inverted."""
    d = gdf.copy()
    for c in (sr, y, la): d[c] = _num(d[c])
    m = (d[sr] > 0) & (d[y] > 0) & (d[la] > 0); f = 0.5 * np.log(d[sr][m] / d[y][m]) + 0.5 * np.log(d[la][m] / d[y][m]); lo, hi = np.percentile(f, [1, 99])
    return minmax(idw_grid(grid, tract, d.geometry.x[m].values, d.geometry.y[m].values, np.clip(np.asarray(f), lo, hi)), tract)

def preservation_from_depth(grid, tract, depth, cover_mask=None, breaks=(100, 300, 600)):
    """Nevada convention (0.74 vs the Nevada layer with a Quaternary cover mask): under cover, depth to basement <=100 m scores 1.0, 100-300 m 0.8,
    300-600 m 0.4, deeper 0.1; ground without cover scores 0.6. Breaks are user parameters (the 600 m limit reflects open-pit depth)."""
    cm = np.ones(depth.shape, bool) if cover_mask is None else cover_mask; b1, b2, b3 = breaks; s = np.full(depth.shape, 0.6)
    s = np.where(cm & (depth <= b1), 1.0, s); s = np.where(cm & (depth > b1) & (depth <= b2), 0.8, s); s = np.where(cm & (depth > b2) & (depth <= b3), 0.4, s); s = np.where(cm & (depth > b3), 0.1, s)
    return np.where(tract & np.isfinite(depth), s, np.nan)

def load_geochron(ar_csv=None, kar_csv=None, samples_csv=None, crs="EPSG:26911", lonlat_crs="EPSG:4326"):
    """USGS geochronology release: join age tables (Ar-Ar, K-Ar) to the sample table (lat/long), keep one age per sample (smallest uncertainty),
    drop -9999 placeholders. Returns a GeoDataFrame with column 'age' (Ma)."""
    import pandas as pd, geopandas as gpd
    parts = []
    for f in (ar_csv, kar_csv):
        if f:
            d = pd.read_csv(f, usecols=lambda c: c in ("sampleid", "age", "ageunc"), low_memory=False); parts.append(d)
    a = pd.concat(parts); a["age"] = pd.to_numeric(a["age"], errors="coerce"); a = a[a["age"] > -9000].dropna(subset=["age"])
    a["ageunc"] = pd.to_numeric(a.get("ageunc"), errors="coerce").fillna(1e9); a = a.sort_values("ageunc").drop_duplicates("sampleid")
    s = pd.read_csv(samples_csv, usecols=["sampleid", "lat", "long"], low_memory=False).dropna(); d = a.merge(s, on="sampleid")
    g = gpd.GeoDataFrame(d, geometry=gpd.points_from_xy(d["long"], d["lat"]), crs=lonlat_crs); return g.to_crs(crs)

def buried_pluton_proximity(grid, tract, grav, mag, intrusive_mask, winsor_q=98):
    """Approximate recipe from the Nevada report: gravity-low plus magnetic-high anomalies calibrated on known intrusive rocks, geometric-mean
    score, converted to proximity. grav and mag are arrays on the grid; intrusive_mask marks mapped intrusive rocks (calibration)."""
    zg = (grav - np.nanmean(grav[intrusive_mask])) / np.nanstd(grav[intrusive_mask]); zm = (mag - np.nanmean(mag[intrusive_mask])) / np.nanstd(mag[intrusive_mask])
    s = np.sqrt(np.clip(-zg, 0, 3) / 3 * np.clip(zm, 0, 3) / 3); s = np.where(np.isfinite(s), s, 0); cand = s >= np.nanpercentile(s[tract], 90)
    return proximity_score(cand & tract, tract, winsor_q)
