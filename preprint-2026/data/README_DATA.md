# Meteorological data provenance

Daily meteorological observations were obtained from Chilean Meteorological Directorate DMC (RE1018) using INIA-operated stations 180030 Socoroma, 180038 Ticnamar and 180046 Chapiquiña, 2020–2024.

**The source and processed daily observations are intentionally not redistributed here pending the applicable DMC/INIA permission and attribution requirements.**

The public repository includes scripts to retrieve RE1018 monthly pages and build the standardized daily CSV. The user must place the derived, quality-controlled file at `data/clima_diario_tres_localidades_qc.csv` to reproduce the simulations. Confirm the terms for your intended use and cite DMC/INIA appropriately.

Eight Ticnamar readings (22–29 February 2020) around −40°C were flagged as invalid; missing data were not automatically interpolated. Benchmark simulations use the first 100 days of 2021.
