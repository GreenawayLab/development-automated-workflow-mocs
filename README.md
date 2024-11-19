
# Development of an Automated Workflow for Metal-Organic Cage Screening

This repository contains the associated code for data curation and characterization used in the study *Development of an Automated Workflow for Screening the Assembly and Host-Guest Behavior of Metal-Organic Cages (MOCs) Towards Accelerated Discovery*.

The experimental data can be found on [Zenodo](https://zenodo.org).

---

## Directory Structure and Scripts

### 1. Mass Spectrometry (MS) Analysis
**Directory**: `ms`  
**Script**: `moc_ms_analyser.py`

- **Description**: Contains the code and resulting data frames from the MS analysis script.
- **Data Path**: The MS data is located in the `/ms` directory from Zenodo.
- **Usage**: Run the script from the parent directory as follows:
  ```bash
  python moc_ms_analyser.py
  ```
- **Outputs**: The script generates identified masses and formulas with the following fields:
  - `Formula`
  - `found_mz`
  - `theoretical_mz`
  - `found_intensity`
  - `theoretical_intensity`
  - `predicted_charge`
  - `mz_difference`
  - `charge`
  - `splitting_error`

---

### 2. NMR Analysis
**Directory**: `nmr`  
**Script**: `moc_nmr_analyser.py`

- **Description**: Provides code and resulting data frames for NMR data analysis.
- **Data Path**: The NMR data can be found at `/nmr/AB-03-002/01` on Zenodo.
- **Usage**: Run the script from the parent directory, specifying the data path:
  ```bash
  python moc_nmr_analyser.py /nmr/AB-03-002/01
  ```

---

### 3. Opentrons Scripts
**Directory**: `opentrons scripts`

- **Description**: Contains protocols for the Opentrons OT-2 robotic system, including:
  - Automated screening
  - Host-guest assay
  - Solvent benchmarking

---

### 4. Computational Modelling
- **Code Reference**: Computational modeling scripts were adapted from the repository: [MOC-automation](https://github.com/PaulaTeeuwen/MOC-automation.git).

---
