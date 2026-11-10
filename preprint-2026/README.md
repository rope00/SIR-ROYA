# Oregano rust — climate-driven agent-based scenario modeling

Research software and simulation-derived outputs associated with the preprint **Climate-driven agent-based scenario modeling of Puccinia menthae epidemics in oregano-growing areas of northern Chile**.

**Status:** public research-code release accompanying a single-author computational preprint by Nicolás Sáez. The author acknowledges Germán Sepúlveda and Mabel Arismendi for motivating field work. The software does not demonstrate that *Puccinia menthae* has been confirmed in Chile and does not generate calibrated field forecasts.

## Reproduction

Python 3.11+ is recommended. Create an environment and run:

```bash
python -m pip install -r requirements.txt
python code/download_re1018.py --help
python code/simulate_v4.py --input data/clima_diario_tres_localidades_qc.csv --output results/reproduced --full
python code/convergence_runner.py
```

The `convergence_runner.py` script reproduces the seven-case, 10,500-run Monte Carlo convergence experiment with a fresh controlled seed stream. A full source-driven reconstruction is conditional on obtaining daily DMC records. All spatial units are model lattice cells; absolute field risk must not be inferred without calibration.

## Data availability and permissions

The public candidate includes aggregate/model-derived simulation results and scripts, not the DMC/INIA daily station observations. See `data/README_DATA.md`.

## Publication and contributors

**Computational preprint author:** Nicolás Sáez. The original project was coordinated by Germán Sepúlveda, and field observations and photographs were gathered by Germán Sepúlveda and Mabel Arismendi. They reviewed the revised single-author manuscript and, according to their written reply provided to the author, approved the single-author publication and their acknowledgment. No field photographs or third-party daily weather records are distributed here.

The publication PDF will be linked after final arXiv submission. No DOI has yet been assigned. For data provenance and attribution, see `data/README_DATA.md`.

## Reuse rights

The source code and model-derived CSV outputs are **publicly inspectable**; figure-generation scripts are provided. No open-source license has been applied at this stage; applicable copyright therefore remains in place. A reuse license will be chosen only after the rights holder confirms the desired terms. GitHub visibility does not by itself grant permission to copy or redistribute source code.
