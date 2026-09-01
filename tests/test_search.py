"""Guards on the filter-subset search.

Two things can go wrong quietly here. The holdout can leak into the search, in
which case the final number is selected-on and inflated; and a subset's K can be
capped at a range boundary while another's is not, in which case the ranking is
comparing them on unequal terms. Both leave correct-looking output.
"""

import numpy as np
import pytest

from dwarfhunt import search as S


def _fake_mags(n_planets=120, n_z=40, labels=("A", "B", "C", "D")):
    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0, 1, n_planets))
    planet = {l: 1.0 + i * 0.4 * t + rng.normal(scale=.02, size=n_planets)
              for i, l in enumerate(labels)}
    galaxies = []
    for off in (0.0, 0.4):
        z = np.linspace(0, 1, n_z)
        galaxies.append({l: 3.0 + off + i * 0.3 * z + rng.normal(scale=.02, size=n_z)
                         for i, l in enumerate(labels)})
    y = np.concatenate([np.zeros(n_planets), np.ones(2 * n_z)])
    return list(labels), planet, galaxies, y


def test_colour_names_are_n_minus_one_and_adjacent():
    assert S.colour_names(["A", "B", "C", "D"]) == ["A - B", "B - C", "C - D"]


def test_subset_matrix_shape_and_rank():
    labels, planet, galaxies, y = _fake_mags()
    X = S.subset_matrix(labels, planet, galaxies)
    assert X.shape == (len(y), len(labels) - 1)
    assert np.linalg.matrix_rank(X) == len(labels) - 1


def test_holdout_split_is_disjoint_and_complete():
    _, _, _, y = _fake_mags()
    pool, hold = S.holdout_split(y, holdout_frac=0.3, random_state=0)
    assert set(pool).isdisjoint(hold)
    assert len(pool) + len(hold) == len(y)
    assert sorted(np.concatenate([pool, hold])) == list(range(len(y)))
    # stratified: both classes present on both sides
    for part in (pool, hold):
        assert set(np.unique(y[part])) == {0.0, 1.0}


def test_search_never_reads_holdout_rows():
    """Poison the holdout rows; a search that touches them cannot stay finite."""
    labels, planet, galaxies, y = _fake_mags()
    pool, hold = S.holdout_split(y, holdout_frac=0.3, random_state=0)

    poisoned = {l: v.copy() for l, v in planet.items()}
    n_planets = len(next(iter(planet.values())))
    for l in poisoned:
        rows = hold[hold < n_planets]
        poisoned[l][rows] = np.nan

    results = S.search_subsets(labels, poisoned, galaxies, y, sizes=[3],
                               k_candidates=[2, 3], reg_covar=1e-2,
                               search_rows=pool, n_seeds=2, n_init=1)
    assert results, "no subsets evaluated"
    for r in results:
        assert np.isfinite(r["mean"]), "NaN leaked in from the holdout rows"


def test_select_k_flags_a_boundary_hit():
    labels, planet, galaxies, y = _fake_mags()
    X = S.subset_matrix(labels, planet, galaxies)
    # A single-value range can only ever return its own boundary.
    _k, _clf, at_edge = S.select_k(X, y, [4], reg_covar=1e-2, n_init=1)
    assert at_edge
    k, _clf, _edge = S.select_k(X, y, [2, 3, 4, 5, 6, 7], reg_covar=1e-2, n_init=1)
    assert 2 <= k <= 7


def test_search_accepts_a_list_y(tmp_path):
    """search_rows is an index array, so a plain list y used to die with
    "only integer scalar arrays can be converted to a scalar index"."""
    labels, planet, galaxies, y = _fake_mags(n_planets=40, n_z=15)
    out = S.search_subsets(labels, planet, galaxies, list(y), sizes=[3],
                           k_candidates=[2, 3], reg_covar=1e-3,
                           search_rows=np.arange(50), n_seeds=2, n_init=1)
    assert out and 0.0 <= out[0]["mean"] <= 1.0


def test_wavelength_order_is_checked_where_wavelengths_exist():
    """colour_names pairs BY POSITION, so an unordered list silently flips the
    sign of some colours. search.py holds bare labels and cannot see that, so
    the guard lives in planets, on the full filter names."""
    from dwarfhunt import planets

    ordered = ("JWST/MIRI.F1065C", "JWST/MIRI.F1140C", "JWST/MIRI.F1550C")
    assert planets.assert_wavelength_ordered(ordered) == list(ordered)

    with pytest.raises(ValueError, match="not in wavelength order"):
        planets.assert_wavelength_ordered(ordered[::-1])


def test_wavelength_order_uses_the_bandpass_not_the_name():
    """WISE.W3 is nominally "12 um" and sorts to the red end by name, but its
    bandpass means 12.8 um, which places it between F1140C and F1550C."""
    from dwarfhunt import planets

    mixed = ("JWST/MIRI.F1065C", "JWST/MIRI.F1140C", "JWST/MIRI.F1550C",
             "WISE/WISE.W3")
    with pytest.raises(ValueError, match="not in wavelength order"):
        planets.assert_wavelength_ordered(mixed)

    fixed = planets.put_filters_in_wavelength_order(mixed)
    assert planets.assert_wavelength_ordered(fixed) == list(fixed)
    assert fixed.index("WISE/WISE.W3") == 2


def test_winners_curse_grows_with_candidate_count():
    a = S.winners_curse(1, 0.02, 10, trials=400)
    b = S.winners_curse(50, 0.02, 10, trials=400)
    c = S.winners_curse(2000, 0.02, 10, trials=400)
    assert abs(a) < 1e-3, "best-of-1 cannot be inflated"
    assert b > a and c > b


def test_holdout_evaluation_scores_only_holdout_rows():
    labels, planet, galaxies, y = _fake_mags()
    pool, hold = S.holdout_split(y, holdout_frac=0.3, random_state=0)
    out = S.evaluate_on_holdout(labels, planet, galaxies, y, pool, hold,
                                k_candidates=[2, 3, 4], reg_covar=1e-2, n_init=1)
    assert 0.0 <= out["holdout_score"] <= 1.0
    assert out["subset"] == tuple(labels)
