# Repository Overview

This repository contains modifications to LAMMPS to implement constant-electrochemical-potential Langevin dynamics and trajectory reweighting, together with the same version of DeePMD-kit used for TP-classical MD.

## Directory Structure

### 1. `deepmd-kit/source`

* Contains the source code of [DeePMD-kit](./deepmd-kit/source/) with custom modifications.

### 2. `lammps-constP-langevin-reweighting`

* Contains the modified source code of [LAMMPS](./lammps-constP-langevin-reweighting/), based on the stable version released on 2 August 2023 with additional updates.

### 3. `examples`

* Provides working [examples](./examples/) based on analytical model, including input files.

