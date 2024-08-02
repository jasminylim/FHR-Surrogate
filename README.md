# FHR Hybrid Surrogate Model for Fluoride-salt-cooled High-temperature Reactor (FHR)

## Decription
The FHR Surrogate Model is a low-fidelity model of an FHR system with a gFHR core and a two-loop heat transportation system. This repo contains the work written in the manuscipt "Hybrid Surrogate Modeling Framework for the Digital Twin of a Fluoride-salt-cooled High-temperature Reactor (FHR)".

## Surrogate Models
* **systemSurrogatePump** : gFHR system surrogate model given the target power 
  * systemSurrogate : gFHR system surrogate with detailed system state
  * systemSurrogatePump : gFHR system surrogate with detailed system state with pump degradation surrogate model

  For usage details, please refer to the [readme](systemSurrogatePump/readme.md).

## Requirements
* python3
* NumPy
* Pandas (for case tutorials)
* Matplotlib (for case tutorials)
* Jupyter Notebooks (for case tutorials)

# Quick Guide for Python Environment Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy
pip install pandas 
pip install matplotlib
pip install jupyter
```