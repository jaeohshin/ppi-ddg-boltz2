"""Machine-specific inputs for the producer-swap test, overridable by environment.

Nothing here ships with the repository: these are the authors' Zenodo release, the
prepared cohort built from it, the locked mmCSM-PPI blind experiment, and a local
Boltz extraction. Set the environment variables to point at your own copies.

  PREFOLD_EXPERIMENTS  directory holding boltz_clean_baseline_20261101/,
                       boltz_error_audit_20261101/ and mmcsm_ridge_blind_20261115/
  PREFOLD_AUTHOR_PKLS  directory holding the authors' skempi_{mut,wt}.pkl
  PREFOLD_LOCAL_PKLS   directory holding <variant>_{mut,wt}.pkl and the join maps
  VARIANT              which local extraction to compare (default: the matched one)
"""
import os
from pathlib import Path

EXP = Path(os.environ.get('PREFOLD_EXPERIMENTS', '/data/work/pythia_ppi/experiments'))
AUTHOR = Path(os.environ.get('PREFOLD_AUTHOR_PKLS', '/data/work/PreFold/processed_data'))
LOCAL = Path(os.environ.get('PREFOLD_LOCAL_PKLS', '/data/work/local_refold_recycle3'))
VARIANT = os.environ.get('VARIANT', 'scaled_v1_r10_tplfull')

BASE = EXP / 'boltz_clean_baseline_20261101'
PREP = BASE / 'prepared_v1/cohort.csv'
DELTA = BASE / 'ridge_full_v1/delta_features.f32'
SKEMPI = EXP / 'boltz_error_audit_20261101/inputs/syntax/data/skempi_v2.csv'
BLIND_COHORT = EXP / 'mmcsm_ridge_blind_20261115/cohort.csv'

WIDTHS = [147456, 147456, 128]
D = sum(WIDTHS)
OFFSETS = [0, WIDTHS[0], WIDTHS[0] + WIDTHS[1]]
N_PREPARED_ROWS = 5729

def join_map():
    return LOCAL / f'join_map_{VARIANT}.json'

def out_dir():
    d = LOCAL / f'run_{VARIANT}'
    d.mkdir(parents=True, exist_ok=True)
    return d

def missing():
    """Which required inputs are absent, for a clear message instead of a traceback."""
    return [str(p) for p in (PREP, DELTA, SKEMPI, BLIND_COHORT,
                             LOCAL / f'{VARIANT}_mut.pkl', LOCAL / f'{VARIANT}_wt.pkl')
            if not p.exists()]
