"""Guards on the cached magnitude tables.

The failure this protects against is silent and expensive: reuse magnitudes that
were computed against a DIFFERENT planet sample and nothing raises. The arrays
are still the right length, the column names still match, and every colour,
every fit and every reported accuracy downstream is quietly wrong. Same shape as
the config-seam bug, one layer up.

So the cache key carries every input that changes the sample, the full key is
written into the file, and it is re-checked on read rather than trusting the
filename hash.
"""

import json

import numpy as np
import pytest

from dwarfhunt import paths
from dwarfhunt.photometry import (_key_digest, _load_entry, _save_entry,
                                  galaxy_magnitudes, galaxy_magnitudes_swire,
                                  planet_magnitudes)

SAMPLE = dict(num_samples=12, radius_range=(0.6, 1.3), distance=10, rng=2)
F1, F2, F3 = ("JWST/MIRI.F1065C", "JWST/MIRI.F1140C", "JWST/MIRI.F1550C")
TEMPLATE = paths.galaxy_template("K15_templates/MIR_library/MIR0.0.txt")
SWIRE_TEMPLATE = paths.galaxy_template("swire-library/Sc_template_norm.sed")
J, KS = "2MASS/2MASS.J", "2MASS/2MASS.Ks"


def test_cache_hit_returns_identical_arrays(tmp_path):
    a = planet_magnitudes("sonora-bobcat", [F1], cache_dir=tmp_path, verbose=False, **SAMPLE)
    b = planet_magnitudes("sonora-bobcat", [F1], cache_dir=tmp_path, verbose=False, **SAMPLE)
    for key in ("teff", "logg", "feh", "radius", "abs_mag_F1065C"):
        assert np.array_equal(a[key], b[key]), key


def test_adding_a_filter_keeps_the_same_objects(tmp_path):
    """The whole point: a filter added later must describe the same sample."""
    a = planet_magnitudes("sonora-bobcat", [F1], cache_dir=tmp_path, verbose=False, **SAMPLE)
    b = planet_magnitudes("sonora-bobcat", [F1, F2], cache_dir=tmp_path, verbose=False, **SAMPLE)

    for key in ("teff", "logg", "feh", "radius"):
        assert np.array_equal(a[key], b[key]), f"sample changed in {key}"
    assert np.array_equal(a["abs_mag_F1065C"], b["abs_mag_F1065C"])
    assert "abs_mag_F1140C" in b


def test_different_sample_parameters_do_not_share_an_entry(tmp_path):
    a = planet_magnitudes("sonora-bobcat", [F1], cache_dir=tmp_path, verbose=False, **SAMPLE)
    other = dict(SAMPLE, num_samples=SAMPLE["num_samples"] + 1)
    b = planet_magnitudes("sonora-bobcat", [F1], cache_dir=tmp_path, verbose=False, **other)

    assert len(a["teff"]) != len(b["teff"])
    assert len(list(tmp_path.glob("planets_*.npz"))) == 2


@pytest.mark.parametrize("changed", [
    {"rng": 99}, {"distance": 20.0}, {"radius_range": (0.5, 1.3)},
])
def test_every_sample_parameter_changes_the_key(changed):
    """A parameter missing from the key is exactly how wrong magnitudes get reused."""
    base = {"version": 1, "kind": "planets", "model": "sonora-bobcat",
            "num_samples": 12, "radius_range": [0.6, 1.3], "distance": 10.0,
            "rng": 2, "deny": None}
    other = dict(base)
    for k, v in changed.items():
        other[k] = list(v) if isinstance(v, tuple) else v
    assert _key_digest(base) != _key_digest(other), changed


def test_equivalent_integer_seeds_reach_one_entry(tmp_path):
    """rng=2 and rng=np.int64(2) are the same sample and must share a cache entry.

    _key_digest serialises with default=str, so a numpy integer used to be
    stored as the string "2" while the requested key still held np.int64(2).
    The two never compared equal, and every read raised the "written for a
    different sample" error above -- about a sample that was in fact identical.
    """
    seeded = {k: v for k, v in SAMPLE.items() if k != "rng"}
    a = planet_magnitudes("sonora-bobcat", [F1], rng=2, cache_dir=tmp_path,
                          verbose=False, **seeded)
    b = planet_magnitudes("sonora-bobcat", [F1], rng=np.int64(2), cache_dir=tmp_path,
                          verbose=False, **seeded)

    assert np.array_equal(a["teff"], b["teff"])
    assert len(list(tmp_path.glob("planets_*.npz"))) == 1, "seeds split the cache"


@pytest.mark.parametrize("bad", [np.random.default_rng(2), None])
def test_a_seed_that_cannot_name_a_sample_is_refused(bad, tmp_path):
    """A Generator stringifies with its memory address, so it hashes differently
    every process: the cache would never hit, a fresh random sample would be
    drawn on every run, and the numbers would move with no signal at all."""
    seeded = {k: v for k, v in SAMPLE.items() if k != "rng"}
    with pytest.raises((TypeError, ValueError), match="integer seed"):
        planet_magnitudes("sonora-bobcat", [F1], rng=bad, cache_dir=tmp_path,
                          verbose=False, **seeded)


def test_key_mismatch_on_read_raises(tmp_path):
    """Filename hashes are a convenience; the stored key is the source of truth."""
    path = tmp_path / "entry.npz"
    _save_entry(path, {"kind": "planets", "num_samples": 12},
                {"teff": np.arange(3.0)})

    assert _load_entry(path, {"kind": "planets", "num_samples": 12}) is not None
    with pytest.raises(ValueError, match="different sample"):
        _load_entry(path, {"kind": "planets", "num_samples": 13})


def test_coverage_guard_fires_even_on_a_cache_hit(tmp_path):
    """A warm cache must not become a way to smuggle an unsupported filter in."""
    z = np.linspace(0.5, 2, 6)
    galaxy_magnitudes(TEMPLATE, [F1], redshifts=z, cache_dir=tmp_path, verbose=False)
    with pytest.raises(ValueError, match="not covered"):
        galaxy_magnitudes(TEMPLATE, ["2MASS/2MASS.Ks"], redshifts=z,
                          cache_dir=tmp_path, verbose=False)


def test_swire_cache_hit_returns_identical_arrays(tmp_path):
    z = np.linspace(0.0, 2.0, 5)
    a = galaxy_magnitudes_swire(SWIRE_TEMPLATE, [J], redshifts=z,
                                cache_dir=tmp_path, verbose=False)
    b = galaxy_magnitudes_swire(SWIRE_TEMPLATE, [J], redshifts=z,
                                cache_dir=tmp_path, verbose=False)
    assert np.array_equal(a["redshift"], b["redshift"])
    assert np.array_equal(a["mag_J"], b["mag_J"])
    assert np.all(np.isfinite(a["mag_J"]))


def test_swire_adding_a_filter_reuses_the_entry(tmp_path):
    z = np.linspace(0.0, 2.0, 5)
    a = galaxy_magnitudes_swire(SWIRE_TEMPLATE, [J], redshifts=z,
                                cache_dir=tmp_path, verbose=False)
    b = galaxy_magnitudes_swire(SWIRE_TEMPLATE, [J, KS], redshifts=z,
                                cache_dir=tmp_path, verbose=False)
    assert np.array_equal(a["mag_J"], b["mag_J"])
    assert "mag_Ks" in b
    assert len(list(tmp_path.glob("galaxy_*.npz"))) == 1


def test_swire_and_k15_entries_do_not_collide(tmp_path):
    """Distinct cache 'kind' keeps the two libraries' entries apart."""
    z = np.linspace(0.5, 2.0, 5)
    galaxy_magnitudes_swire(SWIRE_TEMPLATE, [F1], redshifts=z,
                            cache_dir=tmp_path, verbose=False)
    galaxy_magnitudes(TEMPLATE, [F1], redshifts=z,
                      cache_dir=tmp_path, verbose=False)
    assert len(list(tmp_path.glob("galaxy_*.npz"))) == 2


def test_swire_coverage_guard_fires_even_on_a_cache_hit(tmp_path):
    z = np.array([10.0, 20.0])
    with pytest.raises(ValueError, match="not covered"):
        galaxy_magnitudes_swire(SWIRE_TEMPLATE, [J], redshifts=z,
                                cache_dir=tmp_path, verbose=False)
