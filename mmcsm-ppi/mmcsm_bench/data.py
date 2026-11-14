"""Loading the three files mmCSM-PPI releases, with the two traps they contain.

    https://biosig.lab.uq.edu.au/mmcsm_ppi/data/train_set
    https://biosig.lab.uq.edu.au/mmcsm_ppi/data/blindtest
    https://biosig.lab.uq.edu.au/mmcsm_ppi/data/blindtest_zemu

Trap 1: in `train_set` the mutation list is unquoted, so a row with three
mutations has three more fields than a row with two. Read from the right.

Trap 2: in `train_set` the `pdb` column is the PDB id for forward entries and
`<pdb id>-<mutations>` for the modelled reverse ones. Counting that column
directly gives 344 "complexes" where there are 130.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

FILES = {
    "train_set": "https://biosig.lab.uq.edu.au/mmcsm_ppi/data/train_set",
    "blindtest": "https://biosig.lab.uq.edu.au/mmcsm_ppi/data/blindtest",
    "blindtest_zemu": "https://biosig.lab.uq.edu.au/mmcsm_ppi/data/blindtest_zemu",
}


@dataclass
class Row:
    pdb: str              # four-character PDB id, lowercase
    mutations: tuple[str, ...]
    ddg: float            # experimental, kcal/mol
    predictions: dict[str, float]

    @property
    def key(self) -> tuple[str, tuple[str, ...]]:
        return self.pdb, tuple(sorted(self.mutations))


def pdb_id(field: str) -> str:
    """`1a22` and `1a22-db319a_eb320a` both name the complex 1a22."""
    return field.split("-")[0].strip().lower()


def _rows(path: Path) -> list[list[str]]:
    with open(path, newline="") as fh:
        return [r for r in csv.reader(fh) if r]


def load_train(path: Path) -> list[Row]:
    out = []
    for r in _rows(path)[1:]:
        if len(r) < 5:
            continue
        out.append(Row(pdb=pdb_id(r[0]),
                       mutations=tuple(x for x in r[1:-3] if x),
                       ddg=float(r[-2]),
                       predictions={"mmCSM-PPI": float(r[-1]), "_type": r[-3]}))
    return out


def load_blind(path: Path) -> list[Row]:
    """595 mutations of four or more point mutations, with three methods' predictions."""
    out = []
    for r in _rows(path)[1:]:
        if len(r) < 6:
            continue
        out.append(Row(pdb=pdb_id(r[0]),
                       mutations=tuple(m for m in r[1].split(",") if m),
                       ddg=float(r[2]),
                       predictions={"FoldX": float(r[3]),
                                    "Discovery Studio": float(r[4]),
                                    "mmCSM-PPI": float(r[5])}))
    return out


def load_zemu(path: Path) -> list[Row]:
    """272 mutations over 24 complexes, against ZEMu."""
    out = []
    for r in _rows(path)[1:]:
        if len(r) < 5:
            continue
        out.append(Row(pdb=pdb_id(r[0]),
                       mutations=tuple(m for m in r[1].split(",") if m),
                       ddg=float(r[2]),
                       predictions={"ZEMu": float(r[3]), "mmCSM-PPI": float(r[4])}))
    return out


def load_predictions(path: Path, name: str) -> dict[tuple[str, tuple[str, ...]], float]:
    """A CSV of your own predictions: columns pdb, mutations, prediction.

    `mutations` is a comma- or semicolon-separated list in the release's own
    notation, matched order-insensitively against the benchmark rows.
    """
    out = {}
    with open(path, newline="") as fh:
        for rec in csv.DictReader(fh):
            muts = tuple(sorted(m.strip() for m in
                                rec["mutations"].replace(";", ",").split(",") if m.strip()))
            out[(pdb_id(rec["pdb"]), muts)] = float(rec["prediction"])
    if not out:
        raise ValueError(f"no rows read from {path} for {name}")
    return out
