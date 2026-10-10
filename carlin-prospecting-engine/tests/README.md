# Tests (need the Terra-X project files; set FA to the project folder)

- `test_builders_vs_nevada.py`: rebuilds Nevada layers from the Nevada raw data and reports rank agreement with the layers used in the Nevada model (`builder_vs_nevada_results.csv`).
- `test_builders_openmineral.py`: the newer builders (alteration, preservation, fertility, ore-fluid footprint, buried plutons, raw USGS geochronology tables) on the OpenMineral raw data, read-only (`builders_openmineral_results.csv`). Also checks that a corrupted fault file is detected.
- `test_end_to_end_nevada.py`: builds all 19 layers from raw Nevada data and runs the zero-shot step. A smoke test of the chain, not a validation (Nevada deposits are in the reference).
- `test_fixed_scale_nevada.py`: rebuilds the twelve fixed-scale layers from Nevada raw data in area-relative and fixed mode and asserts they are identical (max difference 6e-8).
- `archetype_cv_results.csv`: within-family spatial-block cross-validation and cross-family transfer for the three archetypes (made by `tools/build_reference_archetypes.py`).
- `yukon_trial/`: the Yukon foreign-district test (scripts to prepare inputs and build layers, the evaluation outputs, `evaluation_conditions.csv`, `evaluation_negative_control.csv`). Needs the OpenMineral Yukon_Prospectivity folder (read-only).
- Tools: `tools/calibrate_nevada_scale.py` writes `reference/nevada_scale.json`; `tools/build_reference_archetypes.py` writes the archetype references.

## Paths
The tests and tools read the project data from the folder named by the `FA` environment variable (default: the current folder), and the Yukon scripts read `YUKON` (default `Yukon_Prospectivity`). Set these before running, for example `FA=/path/to/project python tests/test_fixed_scale_nevada.py`.

## Final Yukon run
`tests/yukon_trial/final_10_layer_run/` holds the scripts and config of the final 10-layer zero-shot run (`build_v3_layers.py` builds the two 88-102 Ma intrusion layers and the whole-Yukon geochemistry layer, `capture10.py` scores the 28 occurrences, `engine_with_tie_fix.py` gives tied scores their average rank). Its result: 16 of 28 occurrences in the top 10% and 21 of 28 in the top 20% (mean percentile 0.845), against 8-layer 6 and 19 (0.79) and the simple baseline 22 and 23 (0.884). Needs the Yukon Geological Survey open data.
