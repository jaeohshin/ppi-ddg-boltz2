#!/bin/bash
# Producer-swap test for one local Boltz extraction.
#   ./run_variant.sh [variant]        default: scaled_v1_r10_tplfull
# Inputs are located by paths.py; override them with PREFOLD_EXPERIMENTS,
# PREFOLD_AUTHOR_PKLS and PREFOLD_LOCAL_PKLS. PY selects the interpreter.
set -eu
export VARIANT="${1:-scaled_v1_r10_tplfull}"
export PYTHONNOUSERSITE=1
PY="${PY:-python3}"
cd "$(dirname "$0")"
echo "=== join: $VARIANT";      $PY build_join.py
echo "=== compare: $VARIANT";   $PY compare_sources.py
echo "=== diagnose: $VARIANT";  $PY diagnose.py
echo "=== blocks: $VARIANT";    $PY block_arms.py
echo "=== surface: $VARIANT";   $PY weight_swap.py
