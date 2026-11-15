"""The six 3SE3 rows, with the affinities actually in SKEMPI 2.0."""
import pytest
from mmcsm_bench.disambiguate import ddg, resolve

# Kd in nM, wild type and mutant, for the two 3SE3 interfaces.
CANDIDATES = {
    "3SE3_B_A": {"NB44E,SB47Q,YB43H": [ddg(1.5e-6, 2.4e-8)],
                 "NB44A,SB47A,YB43A": [ddg(1.5e-6, 4.1e-7), ddg(1.5e-6, 7.5e-7)]},
    "3SE3_B_C": {"NB44E,SB47Q,YB43H": [ddg(3.4e-7, 4.7e-7)],
                 "NB44A,SB47A,YB43A": [ddg(3.4e-7, 4.6e-7)]},
}


def candidates_for(mutation):
    return {k: v[mutation] for k, v in CANDIDATES.items()}


def test_sign_convention_is_the_negative_of_skempis():
    # A mutation that weakens binding (Kd up) is a positive SKEMPI ddG and a
    # negative mmCSM one.  Getting this backwards matches nothing at all.
    assert ddg(1e-7, 1e-6) < 0
    assert ddg(1e-6, 1e-7) > 0


@pytest.mark.parametrize("label,direction,expected", [
    (ddg(3.4e-7, 4.7e-7), "forward", "3SE3_B_C"),   # row 1038
    (-ddg(3.4e-7, 4.7e-7), "reverse", "3SE3_B_C"),  # row 1039, same pair reversed
    (ddg(1.5e-6, 2.4e-8), "forward", "3SE3_B_A"),   # row 1046
])
def test_each_row_resolves_to_one_interface(label, direction, expected):
    assert resolve(candidates_for("NB44E,SB47Q,YB43H"), label, direction) == expected


def test_a_candidate_standing_for_several_rows_still_resolves():
    # 3SE3_B_A carries two measurements of the alanine triple; either one is
    # enough to identify the interface.
    c = candidates_for("NB44A,SB47A,YB43A")
    assert resolve(c, ddg(1.5e-6, 7.5e-7), "forward") == "3SE3_B_A"
    assert resolve(c, ddg(1.5e-6, 4.1e-7), "forward") == "3SE3_B_A"
    assert resolve(c, ddg(3.4e-7, 4.6e-7), "forward") == "3SE3_B_C"


def test_no_match_is_dropped_not_guessed():
    c = candidates_for("NB44E,SB47Q,YB43H")
    assert resolve(c, 99.0, "forward") is None
    # Close is not good enough: the rule requires an exact match.
    assert resolve(c, ddg(3.4e-7, 4.7e-7) + 0.01, "forward") is None


def test_an_ambiguous_label_is_dropped():
    assert resolve({"a": [1.0], "b": [1.0]}, 1.0, "forward") is None


def test_direction_is_validated():
    with pytest.raises(ValueError):
        resolve({"a": [1.0]}, 1.0, "backward")
