import sys, json, os, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","scripts")); import build_layers
H=os.path.abspath("."); I=f"{H}/inputs"; Y=os.environ.get("YUKON","Yukon_Prospectivity")
W={"thrust known":0.5,"thrust inferred":0.35,"known fault":1.0,"inferred fault":0.7,"concealed fault":0.4}
F=f"{I}/faults_classed_scale.gpkg"; harm={"scale_field":"SCALE","common_scale":250}
def cfg(harmonise=True, tract=f"{H}/tract_east_of_tintina.gpkg", basis="North American margin strata east of the Tintina fault (approximated by splitting the Yukon outline along the mapped Tintina fault)"):
    h = harm if harmonise else {}
    return {"crs":"EPSG:3578","cell":1000,"scaling":"fixed","tract":{"file":tract,"bounds":[31174.59,611111.25,975174.59,1707111.25]},"tract_is_permissive":True,"tract_basis":basis,"layers":{
     "arch_fault_density":{"file":F,"weight_field":"F_CLASS","weights":W,**h},
     "arch_fault_intersection_density":{"file":F,"field":"F_CLASS","a_regex":"^(known|inferred|concealed) fault","b_regex":"thrust",**h},
     "arch_fluid_pathway_connectivity":{"file":F,"radius_km":5,**h},
     "arch_structural_complexity":{"file":F,**h},
     "arch_lithologic_contact_proximity":{"file":f"{I}/reactive_cap_contacts.gpkg","field":"SEL","regex":"x"},
     "trap_reactive_carbonate":{"file":f"{I}/trap_reactive_carbonate.gpkg","field":"SEL","regex":"x"},
     "trap_impermeable_cap":{"file":f"{I}/trap_impermeable_cap.gpkg","field":"SEL","regex":"x"},
     "trap_fluid_chemistry":{"file":f"{Y}/Regional_Geochemical_Surveys_RGS_INAA_250k.shp","pathfinders":["AS_INA_PPM","SB_INA_PPM","AU_INA_PPB"],"immobiles":["TH_INA_PPM","HF_INA_PPM","SC_INA_PPM"]},
     "eng_eocene_belt_density":{"file":f"{I}/window_igneous_points.gpkg","age_field":"age","age_min":68,"age_max":80},
     "src_eocene_intrusive_proximity":{"file":f"{I}/window_igneous_points.gpkg","age_field":"age","age_min":68,"age_max":80},
     "eng_volcanic_intrusive_proximity":{"file":f"{I}/eng_volcanic_intrusive_proximity.gpkg","field":"SEL","regex":"x"},
     "src_jurassic_cretaceous_intrusive_proximity":{"file":f"{I}/src_jurassic_cretaceous_intrusive_proximity.gpkg","field":"SEL","regex":"x"},
     "src_basement_proximity":{"file":f"{I}/src_basement_proximity.gpkg","field":"SEL","regex":"x"}}}
if __name__=="__main__":
    which=sys.argv[1]; subset=sys.argv[2].split(",") if len(sys.argv)>2 else None
    c=cfg(harmonise=(which!="noharm"))
    if subset: c["layers"]={k:v for k,v in c["layers"].items() if k in subset}
    log=build_layers.build(c,f"{H}/layers_{which}"); print(json.dumps({k:(v.get("result"),v.get("scaling","")[:5],v.get("share_at_zero"),v.get("share_at_one"),v.get("mapping_scale_control")) for k,v in log.items() if k!="_meta"},indent=0)); print(log["_meta"])
