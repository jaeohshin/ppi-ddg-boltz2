"""The two parsing traps, pinned."""
from mmcsm_bench.data import pdb_id


def test_pdb_id_strips_the_modelled_mutation_suffix():
    assert pdb_id("1a22") == "1a22"
    assert pdb_id("1A22") == "1a22"
    # the reverse entries carry their mutations in the pdb column
    assert pdb_id("1a22-db319a_eb320a_kb321a") == "1a22"


def test_unquoted_mutation_lists_are_read_from_the_right(tmp_path):
    from mmcsm_bench.data import load_train
    p = tmp_path / "train.csv"
    p.write_text("pdb,mutations,type,ddg,prediction_mmCSM-PPI\n"
                 "1tm1,RI67A,RI65A,forward,-5.05192,-2.8\n"
                 "1xyz,AA1B,BB2C,CC3D,reverse,1.5,0.9\n")
    rows = load_train(p)
    assert rows[0].pdb == "1tm1" and rows[0].mutations == ("RI67A", "RI65A")
    assert rows[0].ddg == -5.05192
    assert rows[1].mutations == ("AA1B", "BB2C", "CC3D") and rows[1].ddg == 1.5
