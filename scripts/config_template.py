"""
config.py — shared configuration for the teleseismic MT + FFI pipeline.

Copy this file to the parent directory of scripts/ as `config.py` and edit paths.
All pipeline scripts do `import config as C` and use C.DATA_DIR / C.RESULTS_DIR:
  data/mtuq/<TAG>/          prepared MTUQ data (SAC *.z/*.r/*.t + weights.dat)
  results/mtuq_*_<TAG>/     inversion outputs
"""
import os

BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
DATA_DIR    = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
