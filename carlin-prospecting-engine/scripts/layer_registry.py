"""Registry of the 19 region-agnostic Carlin evidence layers (Nevada-derived).
Removed from the original 23: palinspastic Eocene proximity (redundant, Nevada-specific), lead-isotope crustal suture and strontium isotope line
(Nevada-specific), and proximity to intrusions of several ages (Nevada-specific, low value).
Layer ids keep the Nevada file names so the frozen reference matches; the names shown to users are in NAME.
Several layers are PROXIES for a geological idea (alteration footprint, ore-fluid footprint, preservation). The user picks the proxy that fits the
local geology; the Nevada reference was built with the proxy noted in the name."""
LAYERS = [
 # id, name shown to users, mineral-system part, category (A direct / B needs translation / P proxy chosen for the local geology), inputs
 ("arch_fault_density","Fault density","Pathway","A","Fault lines"),
 ("arch_fault_intersection_density","Fault intersection density","Pathway","A","Fault lines split into high-angle conduit faults and thrusts"),
 ("arch_fluid_pathway_connectivity","Fluid pathway connectivity (fault network)","Pathway","A","Fault lines"),
 ("arch_structural_complexity","Structural complexity (fault orientation disorder)","Pathway","A","Fault lines"),
 ("arch_lithologic_contact_proximity","Proximity to selected lithologic contacts","Pathway","A","Geologic map polygons and the contacts that matter locally"),
 ("trap_reactive_carbonate","Reactive host rock (Nevada: ferroan carbonate)","Trap","P","Geologic map and the units the user judges reactive"),
 ("trap_impermeable_cap","Impermeable cap rock","Trap","P","Geologic map and the units that seal"),
 ("trap_argillic_k_depletion","Hydrothermal alteration footprint (Nevada proxy: potassium depletion)","Trap","P","Radiometric potassium (with lithology) or another alteration proxy"),
 ("trap_fluid_chemistry","Ore-fluid geochemical footprint (Nevada proxy: Au/Ag and base metals)","Trap","P","Stream, soil or rock geochemistry"),
 ("trap_preservation","Preservation (cover thickness and erosion level)","Preservation","P","Depth to basement or cover thickness, and a cover mask"),
 ("eng_eocene_belt_density","Intrusive belt density in the mineralization age window","Energy","B","Dated igneous rocks (age, location)"),
 ("eng_volcanic_intrusive_proximity","Proximity to volcanic and intrusive rocks","Energy","B","Geologic map igneous units"),
 ("eng_buried_pluton_proximity","Proximity to buried plutons","Energy","B","Gravity and magnetic data, plus mapped intrusions for calibration"),
 ("eng_mantle_upwelling","Mantle upwelling (deep gravity)","Energy","B","Gravity grid"),
 ("src_eocene_intrusive_proximity","Proximity to intrusions in the mineralization age window","Source","B","Dated intrusions (age, location)"),
 ("src_jurassic_cretaceous_intrusive_proximity","Proximity to older (Jurassic-Cretaceous) intrusions","Source","B","Geologic map intrusive units of the older age"),
 ("src_basement_proximity","Distance from basement outcrops (higher = farther)","Source","B","Geologic map with basement units"),
 ("src_low_magmatic_fertility","Magmatic fertility (Sr/Y and La/Y index)","Source","B","Geochemistry with Sr, Y and La"),
 ("src_crustal_gradient_fav","Crustal gradient favorability (flat gravity)","Source","B","Gravity grid"),
]
IDS=[l[0] for l in LAYERS]
NAME={l[0]:l[1] for l in LAYERS}
PART={l[0]:l[2] for l in LAYERS}
CATEGORY={l[0]:l[3] for l in LAYERS}
INPUTS={l[0]:l[4] for l in LAYERS}
# Spearman agreement between the generic builder and the layer used in the Nevada model (tests/). None = not measured.
AGREEMENT={
 "arch_fault_density":0.98,"arch_fault_intersection_density":0.998,"arch_fluid_pathway_connectivity":0.56,"arch_structural_complexity":0.99,
 "arch_lithologic_contact_proximity":0.16,"trap_reactive_carbonate":0.70,"trap_impermeable_cap":0.96,"trap_argillic_k_depletion":0.66,
 "trap_fluid_chemistry":0.29,"trap_preservation":0.74,"eng_eocene_belt_density":0.995,"eng_volcanic_intrusive_proximity":0.96,
 "eng_buried_pluton_proximity":-0.37,"eng_mantle_upwelling":0.998,"src_eocene_intrusive_proximity":0.80,"src_jurassic_cretaceous_intrusive_proximity":1.0,
 "src_basement_proximity":1.0,"src_low_magmatic_fertility":0.29,"src_crustal_gradient_fav":0.89}
# Layers that are physical (km, km/km2, fabric disorder) and are scaled with fixed Nevada-calibrated cut-offs, so a given score means the same thing in every
# region (reference/nevada_scale.json). All other layers are rescaled inside each study area and are NOT comparable across regions.
FIXED_SCALE = {"arch_fault_density", "arch_fault_intersection_density", "arch_fluid_pathway_connectivity", "arch_structural_complexity", "arch_lithologic_contact_proximity",
               "trap_reactive_carbonate", "trap_impermeable_cap", "eng_eocene_belt_density", "src_eocene_intrusive_proximity", "eng_volcanic_intrusive_proximity",
               "src_jurassic_cretaceous_intrusive_proximity", "src_basement_proximity"}
# Diagnostic evidence = what marks a Carlin-type system specifically (ore-fluid chemistry, alteration). The rest (structure, host rock, intrusions) exists in most
# sedimentary terranes and does not by itself say 'Carlin'.
DIAGNOSTIC = {"trap_fluid_chemistry", "trap_argillic_k_depletion"}
def status(layer_id):
    a = AGREEMENT.get(layer_id)
    if a is None: return "not tested"
    return "matches the Nevada layer" if a >= 0.9 else "approximate" if a >= 0.5 else "substitute (does not reproduce the Nevada layer)"
