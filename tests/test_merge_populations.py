"""Guards on merge_planet_populations.

The failure this exists to prevent is the same shape as the AGN-fraction zip
bug in pop-separation.ipynb's galaxy cell: two populations built from
different filter_names (or a partially warm photometry cache) don't raise on
their own -- np.concatenate on a shorter list of columns just silently keys
the merged dict by whichever tag was checked, and a filter present in one
population but not the other disappears with nothing to notice. So mismatched
columns must raise, naming the tag and what's missing, not merge on the
intersection.
"""

import numpy as np
import pytest

from dwarfhunt.planets import merge_planet_populations


def _population(n, teff_start):
    return {
        "teff": np.arange(teff_start, teff_start + n, dtype=float),
        "logg": np.full(n, 4.5),
        "abs_mag_J": np.linspace(10, 11, n),
    }


def test_merges_in_insertion_order():
    t = _population(3, 575.0)
    y = _population(2, 275.0)

    merged = merge_planet_populations({"sonora-elfowl-t": t, "sonora-elfowl-y": y})

    assert np.array_equal(merged["teff"], np.concatenate([t["teff"], y["teff"]]))
    assert np.array_equal(merged["abs_mag_J"],
                          np.concatenate([t["abs_mag_J"], y["abs_mag_J"]]))


def test_adds_a_source_column_naming_the_originating_tag():
    t = _population(3, 575.0)
    y = _population(2, 275.0)

    merged = merge_planet_populations({"sonora-elfowl-t": t, "sonora-elfowl-y": y})

    assert list(merged["source_model"]) == ["sonora-elfowl-t"] * 3 + ["sonora-elfowl-y"] * 2


def test_source_key_is_configurable():
    t = _population(2, 575.0)
    merged = merge_planet_populations({"sonora-elfowl-t": t}, source_key="family")
    assert "family" in merged and "source_model" not in merged


def test_mismatched_columns_raise_and_name_what_is_missing():
    t = _population(3, 575.0)
    y = _population(2, 275.0)
    del y["abs_mag_J"]

    with pytest.raises(ValueError) as exc:
        merge_planet_populations({"sonora-elfowl-t": t, "sonora-elfowl-y": y})

    message = str(exc.value)
    assert "sonora-elfowl-y" in message
    assert "abs_mag_J" in message


def test_a_source_key_collision_is_rejected():
    t = _population(2, 575.0)
    t["source_model"] = np.array(["already here"] * 2, dtype=object)

    with pytest.raises(ValueError, match="already a column"):
        merge_planet_populations({"sonora-elfowl-t": t})


def test_empty_populations_dict_is_rejected():
    with pytest.raises(ValueError, match="at least one population"):
        merge_planet_populations({})


def test_a_single_population_still_works():
    t = _population(4, 575.0)
    merged = merge_planet_populations({"sonora-elfowl-t": t})

    assert np.array_equal(merged["teff"], t["teff"])
    assert list(merged["source_model"]) == ["sonora-elfowl-t"] * 4
