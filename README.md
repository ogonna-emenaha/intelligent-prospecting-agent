# Terra-X Intelligent Prospecting Agent

**A proprietary, AI-native prospectivity engine developed in-house by Team Terra-X.** It transforms diverse exploration data into geological evidence layers, learns mineral-system signatures from known districts, and maps prospective ground in new regions.

This first version targets **Carlin-type gold** and is anchored on a Nevada reference of 84 deposits.

## What it does

1. **Asks the geology.** You give a study boundary (a permissive tract) and a descriptive model or literature. You approve every choice.
2. **Builds predefined evidence layers.** Up to 19 layers across the mineral-system parts (Source, Energy, Pathway, Trap, Preservation), like an evidence-layer workflow but automated. Each layer states whether it matches Nevada or is a proxy.
3. **Transfers the Nevada archetype.** A positive-unlabeled random forest trained on 40 Nevada magmatic-family deposits scores the new ground (zero-shot, no local deposits needed). Fixed physical scaling keeps values comparable between regions.
4. **Tests itself.** If you have known occurrences, it reports capture against chance and against a simple baseline map.

It reports **prospectivity or similarity, never deposit probability**, and it makes no economic claims.

## Proven on real ground

The agent was built on a Nevada reference of 84 Carlin-type deposits and then applied zero-shot to the Yukon, a different district it had never seen and where no local deposit was used for training. It ran end to end: from raw geological, structural and geochemical files to evidence layers, to a prospectivity map for the whole territory. The method notes in `carlin-prospecting-engine/references/method_notes.md` document the tests.

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

## Scope and roadmap

- **Built for transfer.** The Nevada archetype carries the pathway and trap knowledge to new districts; source and energy layers are built from each district's own geology, using the ore age and the intrusion age that the geologist supplies.
- **Transparent by design.** Every run lists which layers match the Nevada layer and which are documented proxies, so a geologist sees exactly what each layer represents.
- **Next.** Additional districts, additional deposit types through reference packs, and a standalone interface.

## Credits

Developed by Team Terra-X for the Frank Arnott Award. Data used in testing: Nevada public geoscience data and the Yukon Geological Survey open data.

## Licence

Copyright (c) 2026 Team Terra-X. All rights reserved. The repository is public for evaluation and review. See the `LICENSE` file for the terms of use.
