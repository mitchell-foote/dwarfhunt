"""Guards on merge_dwarf_populations.

The failure this exists to prevent is the same shape as the AGN-fraction zip
bug in pop-separation.ipynb's galaxy cell: two populations built from
different filter_names (or a partially warm photometry cache) don't raise on
their own -- np.concatenate on a shorter list of columns just silently keys
the merged dict by whichever tag was checked, and a filter present in one
population but not the other disappears with nothing to notice. So mismatched
PHOTOMETRY columns must raise, naming the tag and what's missing, not merge on
the intersection.

Parameter columns are the other case and get the opposite treatment. Grids
genuinely differ in their axes -- Elf Owl carries c_o_ratio and log_kzz,
Diamondback carries fsed -- and that is a fact about the grids, not a symptom
of a bad call. Those are unioned and NaN-padded so a multi-family population
keeps every axis attached to the rows it describes. The tests below pin both
halves of that split: a missing abs_mag_ still raises, a missing parameter
does not.
"""

import numpy as np
import pytest

from dwarfhunt.dwarfs import merge_dwarf_populations


def _population(n, teff_start):
    return {
        "teff": np.arange(teff_start, teff_start + n, dtype=float),
        "logg": np.full(n, 4.5),
        "abs_mag_J": np.linspace(10, 11, n),
    }


def test_merges_in_insertion_order():
    t = _population(3, 575.0)
    y = _population(2, 275.0)

    merged = merge_dwarf_populations({"sonora-elfowl-t": t, "sonora-elfowl-y": y})

    assert np.array_equal(merged["teff"], np.concatenate([t["teff"], y["teff"]]))
    assert np.array_equal(merged["abs_mag_J"],
                          np.concatenate([t["abs_mag_J"], y["abs_mag_J"]]))


def test_adds_a_source_column_naming_the_originating_tag():
    t = _population(3, 575.0)
    y = _population(2, 275.0)

    merged = merge_dwarf_populations({"sonora-elfowl-t": t, "sonora-elfowl-y": y})

    assert list(merged["source_model"]) == ["sonora-elfowl-t"] * 3 + ["sonora-elfowl-y"] * 2


def test_source_key_is_configurable():
    t = _population(2, 575.0)
    merged = merge_dwarf_populations({"sonora-elfowl-t": t}, source_key="family")
    assert "family" in merged and "source_model" not in merged


def test_mismatched_columns_raise_and_name_what_is_missing():
    t = _population(3, 575.0)
    y = _population(2, 275.0)
    del y["abs_mag_J"]

    with pytest.raises(ValueError) as exc:
        merge_dwarf_populations({"sonora-elfowl-t": t, "sonora-elfowl-y": y})

    message = str(exc.value)
    assert "sonora-elfowl-y" in message
    assert "abs_mag_J" in message


def test_a_source_key_collision_is_rejected():
    t = _population(2, 575.0)
    t["source_model"] = np.array(["already here"] * 2, dtype=object)

    with pytest.raises(ValueError, match="already a column"):
        merge_dwarf_populations({"sonora-elfowl-t": t})


def test_empty_populations_dict_is_rejected():
    with pytest.raises(ValueError, match="at least one population"):
        merge_dwarf_populations({})


def test_a_single_population_still_works():
    t = _population(4, 575.0)
    merged = merge_dwarf_populations({"sonora-elfowl-t": t})

    assert np.array_equal(merged["teff"], t["teff"])
    assert list(merged["source_model"]) == ["sonora-elfowl-t"] * 4


def _elfowl(n, teff_start):
    """A 5-parameter population: the two axes Diamondback does not have."""
    pop = _population(n, teff_start)
    pop["c_o_ratio"] = np.full(n, 0.5)
    pop["log_kzz"] = np.full(n, 7.0)
    return pop


def _diamondback(n, teff_start):
    """A 4-parameter population: the one axis Elf Owl does not have."""
    pop = _population(n, teff_start)
    pop["fsed"] = np.full(n, 2.0)
    return pop


def test_differing_parameter_axes_are_unioned_not_dropped():
    t = _elfowl(3, 575.0)
    d = _diamondback(2, 900.0)

    merged = merge_dwarf_populations({"sonora-elfowl-t": t,
                                       "sonora-diamondback-highres": d})

    # Every axis survives -- this is the whole point of sampling two grids.
    assert {"c_o_ratio", "log_kzz", "fsed"} <= set(merged)
    assert len(merged["teff"]) == 5


def test_a_padded_axis_is_nan_exactly_on_the_grid_that_lacks_it():
    t = _elfowl(3, 575.0)
    d = _diamondback(2, 900.0)

    merged = merge_dwarf_populations({"sonora-elfowl-t": t,
                                       "sonora-diamondback-highres": d})

    # fsed: real for the two Diamondback rows, NaN for the three Elf Owl rows.
    assert np.array_equal(np.isnan(merged["fsed"]),
                          np.array([True, True, True, False, False]))
    assert np.array_equal(merged["fsed"][3:], d["fsed"])

    # log_kzz: the mirror image.
    assert np.array_equal(np.isnan(merged["log_kzz"]),
                          np.array([False, False, False, True, True]))
    assert np.array_equal(merged["log_kzz"][:3], t["log_kzz"])


def test_padding_is_reported_rather_than_silent(capsys):
    t = _elfowl(3, 575.0)
    d = _diamondback(2, 900.0)

    merge_dwarf_populations({"sonora-elfowl-t": t,
                              "sonora-diamondback-highres": d})

    out = capsys.readouterr().out
    assert "fsed" in out and "sonora-elfowl-t" in out
    assert "log_kzz" in out and "sonora-diamondback-highres" in out


def test_photometry_still_raises_even_when_parameters_legitimately_differ():
    """The relaxation must not become a way to smuggle a filter mismatch past."""
    t = _elfowl(3, 575.0)
    d = _diamondback(2, 900.0)
    t["abs_mag_H"] = np.linspace(9, 10, 3)   # only one side got the H filter

    with pytest.raises(ValueError) as exc:
        merge_dwarf_populations({"sonora-elfowl-t": t,
                                  "sonora-diamondback-highres": d})

    message = str(exc.value)
    assert "abs_mag_H" in message
    assert "sonora-diamondback-highres" in message


def test_column_order_is_first_seen_not_hash_order():
    """Iterating the key union as a set orders columns by PYTHONHASHSEED."""
    t = _elfowl(3, 575.0)
    d = _diamondback(2, 900.0)

    merged = merge_dwarf_populations({"sonora-elfowl-t": t,
                                       "sonora-diamondback-highres": d})

    expected = list(t) + [c for c in d if c not in t] + ["source_model"]
    assert list(merged) == expected
