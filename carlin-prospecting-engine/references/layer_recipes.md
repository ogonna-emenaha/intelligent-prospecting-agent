# Layer recipes and build status (19 layers)

Nevada convention: 1 km cells, every layer rescaled 0-1 within the study area, 1 = most favourable. Agreement = Spearman correlation of the layer rebuilt by the generic builder from raw Nevada data against the layer used in the Nevada model (`tests/`). 0.9 or more = matches; 0.5 to 0.9 = approximate; below 0.5 = substitute (the builder is a different definition, so the user must accept it as a proxy).

Layers marked P are proxies for a geological idea. The user picks the proxy that fits the local geology; the Nevada reference used the proxy named.

| Layer | Part | Recipe | Agreement | Status |
|---|---|---|---|---|
| Fault density | Pathway | Weighted line density, 10 km radius (known 1.0, inferred 0.7, concealed 0.4, thrust known 0.5, thrust inferred 0.35) | 0.98 | matches |
| Fault intersection density | Pathway | Distance to nearest intersection of high-angle conduit faults with thrusts, inverted (1 - d/dmax). All-against-all intersections are not the same layer | 0.998 | matches |
| Fluid pathway connectivity (fault network) | Pathway | Total length of the connected fault network each fault belongs to, averaged within 5 km, log-scaled. The Nevada definition was not recovered | 0.56 | approximate |
| Structural complexity | Pathway | 1 minus mean resultant length of fault strike, 10 km radius | 0.99 | matches |
| Proximity to selected lithologic contacts | Pathway | Nevada used a selected subset of contacts; all unit boundaries match only 0.16. Ask which contacts matter | 0.16 | substitute |
| Reactive host rock (P) | Trap | Distance to the units the user names as reactive (Nevada: ferroan carbonate, optimised threshold) | 0.70 | approximate |
| Impermeable cap rock (P) | Trap | Distance to the units the user names as sealing | 0.96 | matches |
| Hydrothermal alteration footprint (P) | Trap | Nevada proxy: potassium lower than expected for the rock type (per-unit baseline, 1-99% winsorised). Without a lithology baseline, a ln K ~ ln Th residual gives only 0.31. Other proxies can be supplied as a raster | 0.66 | approximate |
| Ore-fluid geochemical footprint (P) | Trap | Nevada proxy: Au/Ag above 10 with absent base metals (needs gold assays). Built here as log of pathfinder elements over immobile elements, IDW k=12 | 0.29 | substitute |
| Preservation (P) | Preservation | Under cover, depth to basement up to 100 m = 1.0, 100-300 m = 0.8, 300-600 m = 0.4, deeper = 0.1; ground without cover = 0.6 | 0.74 | approximate |
| Intrusive belt density in the mineralization age window | Energy | Quartic kernel density of dated igneous samples in the user's age window. The Nevada report says 15 km but the delivered raster needs about 140 km, so 140 km is the default | 0.995 (0.974 from the raw USGS tables) | matches |
| Proximity to volcanic and intrusive rocks | Energy | Igneous units, 1 - min(d, P98)/P98 | 0.96 | matches |
| Proximity to buried plutons | Energy | Gravity-low plus magnetic-high anomalies calibrated on mapped intrusions. The rebuilt version does not reproduce the Nevada layer | -0.37 | substitute (experimental) |
| Mantle upwelling | Energy | Gravity upward-continued 8 km | 0.998 | matches |
| Proximity to intrusions in the mineralization age window | Source | Distance to dated intrusions in the window, inverted | 0.80 to 0.86 | approximate |
| Proximity to older (Jurassic-Cretaceous) intrusions | Source | Older intrusive units by map code, distance inverted | 1.00 | matches |
| Distance from basement outcrops | Source | Oldest basement units. The Nevada layer is distance FROM outcrop (higher = farther), despite being called "proximity" | 0.99999 | matches |
| Magmatic fertility | Source | 0.5 ln(Sr/Y) + 0.5 ln(La/Y), 1-99% winsorised, IDW to the grid. High = fertile. In Nevada, deposits sit where this index is low (AUC 0.29) and the model learns that; the layer is not inverted | 0.29 | substitute (different sample set) |
| Crustal gradient favorability | Source | 1 minus rescaled log total horizontal gradient of gravity (flat crust favoured) | 0.89 | approximate |

Removed from the original 23 layers: proximity to Eocene intrusions with palinspastic restoration (redundant with the intrusion-window proximity, r = 0.90, and Nevada-specific), the lead-isotope crustal suture and the strontium isotope line (Nevada-specific), and proximity to intrusions of several ages (Nevada-specific, low value). The layer previously called "low magmatic fertility" is now "magmatic fertility".
