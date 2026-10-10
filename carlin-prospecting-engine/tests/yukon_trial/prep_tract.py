import geopandas as gpd, numpy as np, shapely, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely.geometry import LineString, Polygon, MultiPolygon
from shapely.ops import split, unary_union
Y=os.environ.get("YUKON","Yukon_Prospectivity"); D=f"{Y}/Carlin_engine_run_2026-10-08/inputs_derived"
f=gpd.read_file(f"{Y}/Faults.shp"); t=f[f.NAME.fillna("").str.contains("intina",case=False)]
sa=gpd.read_file(f"{D}/study_area.gpkg").geometry.iloc[0]; v=gpd.read_file(f"{D}/validation_carlin_ygs.gpkg")
pts=shapely.get_coordinates(t.geometry.values); print(len(pts), t.total_bounds)
c=pts.mean(0); u,s,vt=np.linalg.svd(pts-c,full_matrices=False); ax=vt[0]; proj=(pts-c)@ax
order=np.argsort(proj); P=pts[order]; pr=proj[order]
bins=np.linspace(pr.min(),pr.max(),25); idx=np.digitize(pr,bins); line=np.array([np.median(P[idx==i],0) for i in range(1,len(bins)) if (idx==i).any()])
d0=(line[0]-line[1]); d0/=np.linalg.norm(d0); d1=(line[-1]-line[-2]); d1/=np.linalg.norm(d1)
line=np.vstack([line[0]+d0*2.0e6,line,line[-1]+d1*2.0e6]); L=LineString(line)
parts=split(sa,L); print("pieces",len(parts.geoms),[round(p.area/1e6) for p in parts.geoms])
vc=v.geometry.unary_union.centroid; east=[p for p in parts.geoms if p.buffer(1).contains(vc) or p.distance(vc)==0]
print("rackla centroid pieces",len(east)); E=east[0]
gpd.GeoDataFrame(geometry=[E],crs=3578).to_file("tract_east_of_tintina.gpkg",driver="GPKG"); print("tract km2",round(E.area/1e6))
fig,axp=plt.subplots(figsize=(7,8)); gpd.GeoSeries([sa],crs=3578).plot(ax=axp,color="#eee",edgecolor="k"); gpd.GeoSeries([E],crs=3578).plot(ax=axp,color="#cfe3ff",edgecolor="b")
t.plot(ax=axp,color="r",lw=.8); v.plot(ax=axp,color="k",markersize=6); axp.plot(*L.xy,"g--",lw=.8); axp.set_xlim(sa.bounds[0],sa.bounds[2]); axp.set_ylim(sa.bounds[1],sa.bounds[3]); plt.savefig("tract_check.png",dpi=90)
