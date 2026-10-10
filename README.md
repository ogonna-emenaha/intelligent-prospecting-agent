# Terra-X Intelligent Prospecting Agent

**A proprietary, AI-native prospectivity engine developed in-house by Team Terra-X.** It transforms diverse exploration data into geological evidence layers, learns mineral-system signatures from known districts, and maps prospective ground in new regions.

This first version targets **Carlin-type gold** and is anchored on a Nevada reference of 84 deposits.

## What it does

1. **Asks the geology.** You give a study boundary (a permissive tract) and a descriptive model or literature. You approve every choice.
2. **Builds predefined evidence layers.** Up to 19 layers across the mineral-system parts (Source, Energy, Pathway, Trap, Preservation), like an evidence-layer workflow but automated. Each layer states whether it matches Nevada or is a proxy.
3. **Transfers the Nevada archetype.** A positive-unlabeled random forest trained on 40 Nevada magmatic-family deposits scores the new ground (zero-shot, no local deposits needed). Fixed physical scaling keeps values comparable between regions.
4. **Tests itself.** If you have known occurrences, it reports capture against chance and against a simple baseline map.

It reports **prospectivity or similarity, never deposit probability**, and it makes no economic claims.

## Tested on

Nevada (reference, 84 deposits) and a zero-shot transfer to the Yukon. The method notes in `carlin-prospecting-engine/references/method_notes.md` describe the tests, the results and their limits, including that a simple baseline map is always reported next to the agent's map.

## Install and use

**As a Claude skill (recommended).** Add `carlin-prospecting-engine.skill` in Claude's Skills settings, then ask: *"Run the Carlin prospecting engine on my study area."* The agent will first ask for your study boundary and a descriptive model.

**Scripts only.** Python 3.9 or later:

```
pip install -r requirements.txt
python carlin-prospecting-engine/scripts/intake.py DATA_FOLDER
python carlin-prospecting-engine/scripts/build_layers.py config.json OUT/layers
python carlin-prospecting-engine/scripts/run_zero_shot.py --layers OUT/layers --out OUT/zero_shot
```

An example configuration is in `carlin-prospecting-engine/references/example_config.json`. A plain-language guide is in `docs/User_Guide.pdf`.

## Repository layout

| Path | What it holds |
|---|---|
| `carlin-prospecting-engine.skill` | Installable skill package |
| `carlin-prospecting-engine/SKILL.md` | The agent's rules and workflow |
| `carlin-prospecting-engine/scripts/` | Layer builders, zero-shot, local-model and evaluation scripts |
| `carlin-prospecting-engine/reference/` | Frozen Nevada reference models, scaling and archetypes |
| `carlin-prospecting-engine/references/` | Descriptive model, layer recipes, method notes |
| `carlin-prospecting-engine/tests/` | Tests and the Yukon foreign-district trial |
| `docs/` | User guide |

## Limitations

- One foreign district so far (Yukon, 28 occurrences, 22 in one belt). A second district is needed before claiming the Nevada archetype adds value over a simple map.
- Some layers are proxies for the Nevada layer (marked in each run's summary). Geochemistry, gravity, radiometrics, fertility and preservation are scaled inside each study area and are not comparable across regions.
- The reference-domain flag is experimental.
- The tests need the project's raw data and are not runnable from this repository alone.

## Credits

Developed by Team Terra-X for the Frank Arnott Award. Data used in testing: Nevada public geoscience data and the Yukon Geological Survey open data.

## Licence

No licence has been chosen yet, so all rights are reserved by default. Add a `LICENSE` file when you decide how others may use the code.
