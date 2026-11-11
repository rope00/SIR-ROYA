# Oregano rust — climate-driven agent-based scenario modeling

Research software and simulation-derived outputs associated with the preprint **Climate-driven agent-based scenario modeling of Puccinia menthae epidemics in oregano-growing areas of northern Chile**.

**Status:** public research-code release accompanying a single-author computational preprint by Nicolás Sáez. The author acknowledges Germán Sepúlveda and Mabel Arismendi for motivating field work. The software does not demonstrate that *Puccinia menthae* has been confirmed in Chile and does not generate calibrated field forecasts.

## Reproduction

Python 3.11+ is recommended. Create an environment and run:

```bash
python -m pip install -r requirements.txt
python code/download_re1018.py --help  # documentation only; does not download
# Obtain DMC daily RE1018 records under provider's usage conditions:
python code/download_re1018.py          # fetches 180 monthly pages; internet required
python code/clean_re1018_v4.py --input clima_diario_tres_localidades.csv --output data/clima_diario_tres_localidades_qc.csv
python code/simulate_v4.py --input data/clima_diario_tres_localidades_qc.csv --output results/reproduced --smoke-test
python code/simulate_v4.py --input data/clima_diario_tres_localidades_qc.csv --output results/reproduced --full
python code/convergence_runner.py      # 10,500 simulations; computationally demanding
```

The `convergence_runner.py` script reproduces the seven-case, 10,500-run Monte Carlo convergence experiment with a fresh controlled seed stream. A full source-driven reconstruction is conditional on obtaining daily DMC records. All spatial units are model lattice cells; absolute field risk must not be inferred without calibration.

## Data availability and permissions

The public candidate includes aggregate/model-derived simulation results and scripts, not the DMC/INIA daily station observations. See `data/README_DATA.md`.

## Publication and contributors

**Computational preprint author:** Nicolás Sáez. The original project was coordinated by Germán Sepúlveda, and field observations and photographs were gathered by Germán Sepúlveda and Mabel Arismendi. They reviewed the revised single-author manuscript and, according to their written reply provided to the author, approved the single-author publication and their acknowledgment. No field photographs or third-party daily weather records are distributed here.

The publication PDF will be linked after final arXiv submission. The curated source and simulation outputs are permanently archived on Zenodo as software version 1.0.0: https://doi.org/10.5281/zenodo.24142228 (all-versions DOI: https://doi.org/10.5281/zenodo.24142225). For data provenance and attribution, see `data/README_DATA.md`.

## Reuse licenses (2026)

Copyright © 2026 Nicolás Sáez. Licenses are **scoped to this subproject**;
other material elsewhere in `SIR-ROYA` is not affected.

- **Original Python code:** MIT License. See [code/LICENSE](code/LICENSE).
- **Original simulation-derived CSV outputs:** CC BY 4.0 (where protectable). See [results/LICENSE](results/LICENSE).
- **Original research figures:** CC BY 4.0. See [figures/LICENSE](figures/LICENSE).
- **Original preprint text and figures:** CC BY 4.0. The same notice is in the manuscript.

Please cite the software and the eventual arXiv preprint in scholarly use.
This grant does **not** include photographs, external publications, or
meteorological records from DMC/INIA; see [DATA_POLICY.md](DATA_POLICY.md).
See [LICENSES.md](LICENSES.md) for the complete scope and the CC BY legal text.
**Cite the deposited version:** Sáez, N. (2026). *Oregano rust climate-driven agent-based model: Python source and Monte Carlo simulation outputs (v1.0.0)* [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24142228
