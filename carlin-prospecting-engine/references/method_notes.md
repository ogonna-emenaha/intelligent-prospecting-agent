# Method notes: Nevada reference and what was tested

## Reference population
- 84 Nevada deposits in the project file. Clustering gave three genetic families (54 Eocene-magmatic, 12 older/buried-magma, 18 magma-poor).
- 40 of the 54 Eocene-magmatic deposits passed the falsification envelope (robust Mahalanobis, 95% chi-square) and form the reference set. 37 of the 40 fall inside the modelled permissive tract.
- The recipe was reproduced exactly from the project scripts (map correlation 0.99 with the stored 23-layer magmatic map).

## Changes made for transfer
- 23 layers reduced to 19: palinspastic Eocene proximity, lead-isotope suture, strontium isotope line and multi-age intrusion proximity removed.
- The original recipe fits PCA on all 84 deposits. Fitted on the 40 magmatic deposits alone, PCA changed the map completely (Spearman 0.25 against 0.80 when fitted on all 84). The engine therefore fits the forest on the layers directly, with no PCA. In Nevada hold-out tests the choice made little difference (median percentile 0.95 to 0.98).
- The 19-layer engine map correlates 0.58 with the original 23-layer map; the original recipe run on the reduced layers correlates 0.77. Both maps are stored in `reference/`.

## Nevada spatial-block hold-out (5 folds by deposit location; 40 deposits)
Held-out deposits are scored against a random sample of tract cells; AUC is the mean percentile.

| Setup | AUC | top 5% | top 10% |
|---|---|---|---|
| 23 original layers | 0.984 | 0.95 | 1.00 |
| 19 agnostic layers | 0.980 | 0.93 | 0.98 |
| 19 layers, no Energy | 0.967 | 0.80 | 0.93 |
| 19 layers, no Pathway | 0.966 | 0.80 | 0.88 |
| 19 layers, no Trap | 0.972 | 0.90 | 0.95 |
| 19 layers, no Source | 0.982 | 0.95 | 1.00 |
| 19 layers, no Preservation | 0.981 | 0.95 | 1.00 |

Layer dropout (12 random subsets per size):

| Layers used | AUC (mean) | AUC (worst subset) | top 5% | top 10% |
|---|---|---|---|---|
| 15 | 0.975 | 0.962 | 0.89 | 0.95 |
| 12 | 0.972 | 0.954 | 0.85 | 0.96 |
| 10 | 0.970 | 0.948 | 0.85 | 0.95 |
| 8 | 0.955 | 0.891 | 0.73 | 0.90 |
| 6 | 0.915 | 0.692 | 0.55 | 0.75 |
| 4 | 0.914 | 0.778 | 0.46 | 0.71 |
| 3 | 0.847 | 0.692 | 0.32 | 0.55 |

Reading: with 10 or more layers the held-out deposits stay about as well ranked as with all layers, so the provisional rule is at least 10 layers and 3 of the 5 mineral-system parts. Losing the Energy or Pathway layers costs the most.

**Limits of this test.** The magmatic family was defined by clustering on these same layers, so held-out deposits are consistent with the family by construction. Scores near the ceiling show internal consistency, not discovery power. A foreign district with held-out deposits is the real test and has not been run.

## Archetype comparison (diagnostic)
When the user has fewer than 19 layers, the zero-shot step refits the Nevada reference on those layers only (no reweighting or tuning to compensate). It then reports how closely that reduced model matches the full 19-layer archetype on Nevada: rank agreement of the two Nevada maps, overlap of their top 10% cells, and the held-out Nevada deposit ranking of each (5 spatial blocks). It does not alter the user's map and does not measure its accuracy. Example, 10 of 19 layers: map agreement 0.92, top-10% overlap 0.83, held-out deposit AUC 0.979 against 0.980 for the full archetype.

## Findings that change how layers are built
- Intrusive belt density: the report says a 15 km kernel, but the delivered raster is reproduced only by a kernel near 140 km (0.995). The layer measures regional magmatic-belt presence.
- "Proximity to basement" is stored as distance from basement outcrops (higher = farther).
- Lithologic contacts: Nevada used a selected subset of contacts (0.16 with all boundaries).
- Fault intersections must be between high-angle conduit faults and thrusts (0.998; all-against-all gives 0.48).
- Preservation in Nevada is five classes by depth under cover (1.0, 0.8, 0.4, 0.1) and 0.6 elsewhere.
- Potassium depletion needs a lithology baseline (0.66); a regional K-Th regression gives 0.31.
- The layer called "low magmatic fertility" stores the fertility index itself (high = fertile). Deposits sit where it is low. It is renamed magmatic fertility.
- A consistency check on the project's OpenMineral folder found that `Faults_NAD83.gpkg` has every coordinate infinite (54,713 features). It is corrupted and must be re-exported; the project's working copy under Nevada Datasets is valid.

## Reference-domain flag
A kNN distance flag was too sensitive: on Nevada layers rebuilt from raw data it flagged 76% of Nevada, the same as for column-shuffled data. A robust Mahalanobis flag is used instead (about 6% of rebuilt Nevada flagged). It catches gross extrapolation only and has not been shown to separate Carlin from non-Carlin ground.

## Three archetypes (Nevada, tests/archetype_cv_results.csv)
Clusters from the original 84-deposit clustering: magmatic (54-member cluster, 40 inside the falsification envelope), magma-poor (18) and older/buried-intrusion (12). Each gets its own positive-unlabeled forest (`reference/reference_model_<archetype>.joblib`; deposits in `reference_archetypes.npz`). Spatial-block cross-validation inside each family: magmatic mean held-out percentile 0.975 (all 19 layers) and 0.952 (without 5 intrusion-dependent layers); magma-poor 0.832 and 0.848; older/buried 0.889 and 0.861. Cross-family (the tested family is never in training): magmatic model on magma-poor deposits 0.737 (0.796 without intrusion layers) and on older/buried 0.595; magma-poor model on magmatic 0.832 and on older/buried 0.385; older/buried model on magmatic 0.871 and on magma-poor 0.528. The families are genuinely different in Nevada, and magma-poor and older/buried do not transfer to each other. The families were clustered on these same layers, so within-family numbers are optimistic. The magma-poor and older/buried families have only 18 and 12 deposits.

## Fixed physical scaling (tests/test_fixed_scale_nevada.py)
Until now each layer was rescaled 0-1 inside the study area, so the same score meant different distances in different regions. Twelve physical layers now use cut-offs calibrated once on the Nevada tract (`reference/nevada_scale.json`; distance caps in metres, density and fabric ranges). Nevada's own layers were scaled with exactly these cut-offs, so rebuilding Nevada in fixed mode reproduces the area-relative layers (largest difference 6e-8) and the stored reference stays valid. Caps are Nevada's 98th-percentile distance (or the farthest cell for fault intersections): basement 255 km, older intrusions 275 km, intrusion window 171 km, intersections 118 km, volcanic/intrusive 66 km, impermeable cap 48 km, reactive host 40 km, contacts 5 km; fault density top 1.04 km per km2. Layers not in this list (ore-fluid geochemistry, gravity-derived, radiometric potassium, fertility, preservation) remain area-relative and are flagged in every summary. Fault density is reported in km of weighted fault per km2 and the intrusive belt density in samples per km2, so neither depends on cell size.

## Mapping-detail control for faults
`generalise_lines` brings faults compiled at different scales to one detail (features finer than the common scale are simplified at 0.5 mm of map scale and dropped if shorter than 2 km). In the Yukon file, 1,027 of 34,862 faults were dropped at the 1:250,000 common scale. This is an approximation; it removes the largest mapping-intensity artefact but not all.

## Foreign-district test: Yukon (tests/yukon_trial/)
Setup: 28 Carlin-type occurrences from the Yukon Geological Survey (MINFILE) in 6 areas (single-linkage 15 km); 22 are in one belt. Study area: the Yukon outline split along the mapped Tintina fault, keeping the north-eastern North American-margin side (295,000 km2), an approximation. Thirteen layers on a 1 km grid with fixed scaling (geology-derived masks recovered from an earlier run because the unit selections were not saved). Negative control: about 1,100 non-Carlin MINFILE occurrences in the tract in three groups. Percentiles are measured against 30,000 random cells.
1. Baseline (reactive host plus fault density, equal rank weights): mean deposit percentile 0.872, 21 of 28 deposits in the top 10%.
2. Zero-shot: magmatic 0.723, magma-poor 0.255, older/buried 0.755. No archetype beat the baseline. The layer sign check shows why: the magmatic archetype relies on intrusion-window layers where Yukon deposits have the opposite sign (district AUC 0.30 to 0.33 against 0.87 to 0.94 in Nevada) and on basement distance (0.35 against 0.75); magma-poor reverses on fault density (0.84 against 0.37), contacts (0.90 against 0.42) and basement distance (0.35 against 0.91).
3. Few-shot, held-out whole areas, same layers in both arms (mean deposit percentile, local only against Nevada plus local): older/buried one area 0.849 against 0.912 (top 10%: 0 against 22 of 28; areas in the top 20%: 2 against 5 of 6), two areas 0.918 against 0.932, five areas 0.957 against 0.956; magma-poor 0.821 against 0.801, 0.907 against 0.902, 0.972 against 0.962; magmatic 0.888 against 0.872, 0.940 against 0.920, 0.973 against 0.947. Only the older/buried archetype shows a gain with one or two local areas, and it vanishes with five.
4. Negative control (AUC of Carlin against other types; 0.5 is no separation): baseline 0.66 against non-Carlin sediment-hosted (MVT, SEDEX, manto, barite), 0.83 against other gold, 0.83 against other metals; magmatic zero-shot 0.61, 0.38, 0.57; older/buried 0.68, 0.66, 0.81; magma-poor 0.31, 0.08, 0.16. The baseline is therefore specific to Carlin-type only weakly against other sediment-hosted deposits.
5. Geochemistry (an ore-fluid layer built from INAA As, Sb, Au over immobile elements) helped: mean percentile 0.607 to 0.681 for magmatic and 0.742 to 0.768 for older/buried. A diagnostic blend (weight 0.5) raised the share of deposits in the top 20% (0.14 to 0.25; 0.43 to 0.61) and lowered the mean percentile (0.68 to 0.63; 0.77 to 0.73), so the default stays at a single joint forest.
Limits: one district; 22 of 28 deposits in one belt, so held-out areas are not independent of the training areas (a single other area predicted the held-out area well for every model); the older/buried archetype was singled out after seeing the result, so its gain needs a second district; the geology masks and the tract are approximations; gravity layers were left out.

## Yukon follow-up: pathway and trap, then source and energy re-specified (final 10-layer run)
Pathway and trap layers alone (8 layers) ranked best zero-shot with the magmatic archetype (mean percentile 0.91). Source and energy were then re-specified from the Yukon geology table: proximity to 88-102 Ma intrusive polygons (Tombstone, Mayo, Tungsten, Tay River and others; fixed 275 km cap) and their belt density (140 km quartic kernel, fixed range 0 to 0.00352956). These carry a weak signal on their own (AUC about 0.62 and 0.63). Yukon mineralization is about 74 Ma with no mapped coeval district-scale intrusion, so these intrusions are a regional proxy, not coeval with ore.
Final run: 10 layers (5 pathway, 3 trap incl. ore-fluid geochemistry, 1 source, 1 energy), whole Yukon, Nevada magmatic archetype, zero-shot. 16 of 28 occurrences in the top 10% and 21 of 28 in the top 20% (mean percentile 0.845). The earlier 8-layer run gave 6 and 19 (0.79); the simple baseline gives 22 and 23 (0.884). The gain over the earlier run comes from a changed layer set (fluid chemistry added, basement distance removed, fixed scaling), not from the two intrusion layers alone. The simple baseline still ranks the occurrences higher; 22 of 28 occurrences lie in one belt, so the result is indicative only.
Consequence for the skill: pathway and trap transfer from Nevada, source and energy are district-specific, so the agent asks for the ore age and the source-intrusion age separately and builds source and energy from the local geology.

## Still to do
1. A second foreign district (preferably a different belt, with more than one independent area) to confirm or reject the older/buried gain.
2. Match the Nevada definitions for contacts, buried plutons, fluid chemistry (Au/Ag), fertility and connectivity, or replace them with documented proxies. The Yukon sign check suggests ore-fluid chemistry, fault density and basement distance need district-specific checking.
3. Make the physical scaling cover the geochemistry, gravity and radiometric layers (needs a standard unit, for example a z-score against a regional background).
4. Weight the diagnostic evidence (As, Sb, Tl, Hg, decarbonatization, silicification) with a proper Nevada layer; the Nevada reference only has an Au/Ag proxy.
