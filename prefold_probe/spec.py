"""What the reader is, and which experiment settled each choice.

Every entry names the decision and the measurement behind it, so that a reader
of the code does not have to take the design on trust.  Dates refer to the
project ledger in the companion working repository.
"""
from __future__ import annotations

SPEC = {
    "target": (
        "y = log10(Kd_mutant in nM) - log10(Kd_wild-type in nM), the label the "
        "published model is trained and scored on. Multiply by 1.364 for kcal/mol. "
        "The upstream code never converts, so published MAE and RMSE are in label "
        "units, not energies."
    ),
    "input": (
        "The mutant-minus-wild-type difference of three released blocks: wsi and "
        "ws (384x384 each, interface-pooled outer products of the input and trunk "
        "single representations) and wz (128, the interface-pooled pair "
        "representation). Nothing else. No sequence model, no structure beyond "
        "what Boltz-2 already pooled."
    ),
    "model": (
        "Kernel ridge, one penalty scale per block, each block's Gram divided by "
        "its width. Closed form: no seed, no epochs, no early stopping."
    ),
    "why_block_penalty": (
        "A single penalty over the concatenation scores about 0.53, because the "
        "two 147,456-wide matrices dominate the 128-wide wz where most of the "
        "signal is. Per-block scaling gives 0.59 to 0.62. This is the whole "
        "difference between a reader that loses to the head and one that beats it "
        "(2026-10-07 b, c). An earlier reading of this - that block penalties "
        "protect weak blocks from strong ones - was withdrawn (2026-10-08 b): the "
        "0.049 is the cost of the single penalty itself."
    ),
    "why_no_selection": (
        "The inner alpha-selection surface is flat near its top (0.006 across a "
        "300-fold range of alpha) while the outer consequence is ten times larger "
        "(0.064), and nested-CV's smaller inner training folds bias selection "
        "toward over-regularisation. Three nominally identical runs chose three "
        "different alpha sets. Quote the reader-minus-head difference, which is "
        "alpha-independent over alpha 1 to 1000, not a single reader score."
    ),
    "excluded": {
        "ESM-2 650M as an extra block": "0.596 vs 0.615 for the reader alone, not significant (2026-10-07 i)",
        "Pythia-PPI predictions as a block": "0.587 vs 0.593, not significant; output averaging 0.591 (2026-10-08 b)",
        "richer inputs from the saved traces": "response-ranked interface edges, distance shells and earlier blocks give at most +0.01 over a reader that reads all three tensors properly (2026-10-07 a)",
        "other Pairformer layers and recycles": "88 snapshots scanned; flat after recycle 3-4; best +0.005 [-0.019, +0.031] over the final one (2026-09-20)",
        "a nonlinear head on the same input": "MLP minus ridge -0.002 [-0.021, +0.017] on the full data (2026-09-21)",
        "head capacity": "0.1M to 58.4M parameters all score below the reader; a finite tuning grid cannot prove the converse, so the head grid was deliberately not run (2026-10-07 c)",
        "combining with Pythia-PPI": "residual correlation 0.765, and a gate that picks the better model per variant is not learnable (accuracy 0.47 against a 0.50 base rate) (2026-09-30, 2026-10-06)",
    },
    "known_defect": (
        "The published loader pairs a mutant with the first wild-type record "
        "sharing its 4-letter PDB id, so 35 variants in 2C5D, 3SE3 and 3SE4 carry "
        "the wrong wild type. Reproduced here rather than fixed, for comparability "
        "with the published numbers."
    ),
}


def describe() -> str:
    lines = ["prefold-dg linear probe", ""]
    for k in ("target", "input", "model", "why_block_penalty", "why_no_selection", "known_defect"):
        lines += [k.replace("_", " ").upper(), "  " + SPEC[k], ""]
    lines.append("MEASURED AND EXCLUDED")
    for k, v in SPEC["excluded"].items():
        lines.append(f"  {k}: {v}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(describe())
