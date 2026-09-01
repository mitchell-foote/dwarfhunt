"""Guards on the component-to-class mapping.

This is the module that replaced the notebook's `1 - pred` flip, and its logic
is the least visible in the package: the fit never sees a label, so everything
that makes it a *classifier* happens afterwards, in a vote that can be wrong
without being loud. Two ways in particular:

- The vote is over class-BALANCED responsibility mass. Drop the balancing and
  every mixed component swings to whichever class has more training points --
  here the galaxies, ~270 against ~100 dwarfs. The accuracy that comes back is
  still a plausible number.
- A component that draws no responsibility at all has an all-zero mass row, and
  `argmax` on an all-zero row returns 0, silently claiming classes[0]. The
  fallback to the majority class is only observable if you look for it.

Neither leaves a mark in the score. The tests below therefore pin the vote
itself, not just the accuracy it produces, and each one is built so that the
wrong implementation gives a *different answer* rather than a noisier one.

No species import: gmm.py depends only on numpy and sklearn.
"""

import numpy as np
import pytest

from dwarfhunt.gmm import (GMMClassifier, balanced_accuracy,
                           component_class_map)


# --- the class-balanced vote -------------------------------------------------

def test_balancing_lets_a_small_class_win_a_mixed_component():
    """The reason component_class_map divides by class size.

    Component 0 draws responsibility from 30 of 100 class-0 points and 15 of 20
    class-1 points. A raw responsibility sum makes that 30 vs 15 and hands the
    component to class 0; balanced it is 0.30 vs 0.75 and it belongs to class 1.
    The two rules disagree, which is what makes this a test.
    """
    n0, n1 = 100, 20
    y = np.array([0] * n0 + [1] * n1)

    resp = np.zeros((n0 + n1, 2))
    resp[:30, 0] = 1.0             # 30 of the 100 class-0 points
    resp[30:n0, 1] = 1.0
    resp[n0:n0 + 15, 0] = 1.0      # 15 of the 20 class-1 points
    resp[n0 + 15:, 1] = 1.0

    comp_to_class, mass, degenerate = component_class_map(resp, y, np.array([0, 1]))

    assert resp[:, 0].sum() > resp[y == 1, 0].sum(), "raw vote would say class 0"
    assert comp_to_class[0] == 1, "balanced vote should say class 1"
    assert mass[0].tolist() == pytest.approx([0.30, 0.75])
    assert not degenerate.any()


def test_mass_is_responsibility_per_training_point_of_that_class():
    y = np.array([0, 0, 0, 0, 1, 1])
    resp = np.zeros((6, 1))
    resp[:4, 0] = 0.5     # 2.0 spread over 4 class-0 points -> 0.5
    resp[4:, 0] = 1.0     # 2.0 spread over 2 class-1 points -> 1.0

    _comp, mass, _deg = component_class_map(resp, y, np.array([0, 1]))
    assert mass[0].tolist() == pytest.approx([0.5, 1.0])


# --- the degenerate-component fallback ---------------------------------------

def test_a_component_with_no_votes_goes_to_the_majority_not_to_class_zero():
    """argmax on an all-zero row returns 0, which is a silent claim on classes[0].

    Class 1 is deliberately the majority here, so the fallback and the argmax
    bug give different answers.
    """
    y = np.array([0, 0, 1, 1, 1, 1])
    resp = np.zeros((6, 2))
    resp[:, 0] = 1.0               # component 1 draws nothing at all

    comp_to_class, _mass, degenerate = component_class_map(resp, y, np.array([0, 1]))

    assert comp_to_class[1] == 1, "should fall back to the majority class"
    assert degenerate.tolist() == [False, True], "and say that it did"


def test_a_component_with_votes_is_not_flagged_degenerate():
    y = np.array([0, 0, 1, 1])
    resp = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0]])
    _comp, _mass, degenerate = component_class_map(resp, y, np.array([0, 1]))
    assert not degenerate.any()


# --- the classifier on top ---------------------------------------------------

def _blobs(n_per=40, centers=((0, 0), (5, 5), (10, 0)), seed=0):
    rng = np.random.default_rng(seed)
    X = np.vstack([rng.normal(c, 0.35, (n_per, 2)) for c in centers])
    y = np.repeat(np.arange(len(centers)), n_per)
    return X, y


def test_the_fit_never_sees_the_labels():
    """The central claim of the module: fit is unsupervised, labels arrive after."""
    X, y = _blobs()
    a = GMMClassifier(n_components=3, n_init=2, random_state=0).fit(X, y)
    b = GMMClassifier(n_components=3, n_init=2, random_state=0).fit(X, (y + 1) % 3)

    assert np.array_equal(a.gm_.means_, b.gm_.means_)
    assert not np.array_equal(a.comp_to_class_, b.comp_to_class_)


def test_predict_handles_more_than_two_classes():
    """`1 - pred` mapped component 2 to -1; this is what replaced it."""
    X, y = _blobs()
    clf = GMMClassifier(n_components=3, n_init=2, random_state=0).fit(X, y)

    pred = clf.predict(X)
    assert set(np.unique(pred)) <= {0, 1, 2}
    assert clf.score(X, y) > 0.95


def test_more_components_than_classes_still_predicts_real_classes():
    """The intended use: K > n_classes, several blobs tiling one population."""
    X, y = _blobs(centers=((0, 0), (2, 2), (8, 0), (10, 2)))
    y = np.array([0] * 80 + [1] * 80)

    clf = GMMClassifier(n_components=6, n_init=2, random_state=0, reg_covar=1e-3).fit(X, y)

    assert clf.n_components > len(clf.classes_)
    assert set(np.unique(clf.predict(X))) <= {0, 1}
    assert clf.score(X, y) > 0.95


def test_predict_proba_rows_sum_to_one_and_cover_predict():
    """predict is a hard argmax over COMPONENTS, predict_proba a sum over classes,
    so the two may legitimately disagree -- but the class predict chose can never
    be one predict_proba gave zero mass."""
    X, y = _blobs()
    clf = GMMClassifier(n_components=5, n_init=2, random_state=0, reg_covar=1e-3).fit(X, y)

    proba = clf.predict_proba(X)
    assert proba.shape == (len(X), len(clf.classes_))
    assert proba.sum(axis=1) == pytest.approx(np.ones(len(X)))

    chosen = np.searchsorted(clf.classes_, clf.predict(X))
    assert (proba[np.arange(len(X)), chosen] > 0).all()


def test_assert_finite_raises_on_a_corrupted_fit():
    """The positive counterpart to _quiet_blas_matmul: the matmul warnings are
    suppressed as backend noise, so something has to affirm the model is sound
    rather than merely quiet."""
    X, y = _blobs()
    clf = GMMClassifier(n_components=3, n_init=2, random_state=0).fit(X, y)

    clf.gm_.means_[0, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        clf.assert_finite()


def test_component_table_has_one_row_per_component():
    X, y = _blobs()
    clf = GMMClassifier(n_components=4, n_init=2, random_state=0, reg_covar=1e-3).fit(X, y)

    rows = clf.component_table(X, y)
    assert len(rows) == 4
    assert {r["component"] for r in rows} == {0, 1, 2, 3}
    assert sum(r["weight"] for r in rows) == pytest.approx(1.0)
    for cls in clf.classes_:
        assert all(f"n_class_{cls}" in r for r in rows), "raw counts missing"


# --- balanced_accuracy -------------------------------------------------------

def test_balanced_accuracy_is_not_plain_accuracy():
    """4 class-0 and 2 class-1: recalls 3/4 and 1/2 average to 0.625, while
    plain accuracy is 4/6."""
    y_true = np.array([0, 0, 0, 0, 1, 1])
    y_pred = np.array([0, 0, 0, 1, 1, 0])

    assert balanced_accuracy(y_true, y_pred) == pytest.approx(0.625)
    assert np.mean(y_true == y_pred) == pytest.approx(4 / 6)


def test_balanced_accuracy_refuses_to_return_nan():
    """np.mean([]) is nan plus two numpy warnings, and a nan in a score table
    reads as a real number that happens to be missing."""
    with pytest.raises(ValueError, match="none of classes"):
        balanced_accuracy([0, 0, 1], [0, 0, 1], classes=[7, 8])


def test_a_class_missing_from_the_test_fold_is_skipped_not_counted_zero():
    y_true = np.array([0, 0, 0])
    y_pred = np.array([0, 0, 1])
    assert balanced_accuracy(y_true, y_pred, classes=[0, 1]) == pytest.approx(2 / 3)
