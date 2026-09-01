"""Exhaustive search over filter subsets, with the winner's curse handled.

The mechanics are cheap: cached magnitudes make every subset pure arithmetic, so
the whole sweep is seconds to minutes. The statistics are the hard part.

Why the holdout is not optional
-------------------------------
Score N candidates, pick the best, and the winner's score is inflated even when
every candidate is equally good -- you are reporting the maximum of N noisy
draws, not the quality of the winner. Simulated with candidates that are all
truly 0.950 and per-split sigma 0.015:

    candidates    best observed    inflation
             1           0.9500       +0.0000
           165           0.9591       +0.0091
          1981           0.9615       +0.0115
          1981 (5 splits only)         +0.0229

That inflation is comparable to, or larger than, the differences a filter-set
study is trying to resolve. So the search score is a RANKING device and must
never be reported as the winner's performance.

The protocol here enforces that:

1. Objects are split once into a search pool and a holdout.
2. Every subset, every seed, and every K selection happens inside the pool.
3. Finalists are refit on the whole pool and scored on the holdout exactly once.
4. The holdout score is the reportable number.

Step 3 is only honest if the holdout stays untouched until the end. Scoring
finalists, looking, then re-searching spends it.
"""

from itertools import combinations

import numpy as np
from sklearn.model_selection import train_test_split

from .gmm import GMMClassifier
from .planets import color_pairs


def colour_names(labels):
    """Adjacent-pair colour names for a filter subset, blue minus red.

    Assumes `labels` is already in wavelength order. This module cannot verify
    that: it works in bare labels ("F1065C"), which carry no wavelength, and it
    deliberately holds no species dependency -- the whole point of the cached
    magnitude tables is that a subset sweep is pure arithmetic. So the check
    lives one level up, where the full filter names still exist: run the filter
    list through planets.put_filters_in_wavelength_order (or assert it with
    planets.assert_wavelength_ordered) before deriving the labels handed here.

    Out of order, nothing raises -- the colours are still computed and still
    named, but some are the negative of what their name says.

    n filters give n-1 independent colours; adjacent pairs are one such basis.
    """
    return [f"{blue} - {red}" for blue, red in zip(labels, labels[1:])]


def subset_matrix(labels, planet_mags, galaxy_mags):
    """Feature matrix for one filter subset, from cached magnitude dicts.

    planet_mags : {label: ndarray} absolute magnitudes, one entry per filter
    galaxy_mags : list of {label: ndarray}, one dict per template

    Rows stack planets first, then each template in order -- the same row
    convention the notebook's n_planets offset relies on.
    """
    labels = list(labels)
    cols = colour_names(labels)
    planet_cols = color_pairs({l: planet_mags[l] for l in labels}, order=labels)
    galaxy_cols = [color_pairs({l: g[l] for l in labels}, order=labels)
                   for g in galaxy_mags]
    return np.column_stack([
        np.concatenate([planet_cols[c]] + [g[c] for g in galaxy_cols])
        for c in cols])


def select_k(X_fit, y_fit, k_candidates, reg_covar, n_init=10, random_state=0):
    """Pick K by BIC on the data it was fit to. Returns (k, clf, at_edge).

    BIC is computed on the fitting fold, never on the fold used to score -- a
    criterion evaluated on the scoring data is selection on the test set.

    `at_edge` flags a K sitting on a boundary of k_candidates, which means BIC's
    optimum lies outside the range and the returned value is a cap rather than a
    choice. Callers should surface it: a search where some subsets are capped
    and others are not is no longer comparing them on equal terms.
    """
    fits = [GMMClassifier(n_components=k, random_state=random_state,
                          reg_covar=reg_covar, n_init=n_init).fit(X_fit, y_fit)
            for k in k_candidates]
    i = int(np.argmin([f.bic(X_fit) for f in fits]))
    k = k_candidates[i]
    return k, fits[i], k in (k_candidates[0], k_candidates[-1])


def score_subset(X, y, *, k_candidates, reg_covar, n_seeds=10, test_size=0.3,
                 n_init=10, random_state=0):
    """Mean balanced accuracy for one subset over `n_seeds` internal splits.

    Every subset sees the same seeds, hence the same partitions, so the ranking
    is paired across subsets rather than each one drawing its own luck.
    """
    scores = np.empty(n_seeds)
    ks = np.empty(n_seeds, dtype=int)
    edges = 0
    for seed in range(n_seeds):
        Xtr, Xte, ytr, yte = train_test_split(
            X, y, test_size=test_size, stratify=y, random_state=seed)
        k, clf, at_edge = select_k(Xtr, ytr, k_candidates, reg_covar,
                                   n_init=n_init, random_state=random_state)
        ks[seed] = k
        edges += bool(at_edge)
        scores[seed] = clf.score(Xte, yte)
    return {"mean": float(scores.mean()),
            "std": float(scores.std(ddof=1)) if n_seeds > 1 else 0.0,
            "scores": scores, "k": ks, "n_at_edge": edges}


def search_subsets(labels, planet_mags, galaxy_mags, y, *, sizes,
                   k_candidates, reg_covar, search_rows=None, n_seeds=10,
                   test_size=0.3, n_init=10, progress=None):
    """Rank every filter subset of the requested sizes.

    search_rows : index array restricting the search to the search pool. Pass
        the pool from holdout_split; leaving it None searches everything, which
        is fine for exploration but leaves no clean data to report on.
    """
    labels = list(labels)
    # asarray, because search_rows is an index array and `y[search_rows]` below
    # raises "only integer scalar arrays can be converted to a scalar index" on
    # a plain list -- a confusing error a long way from the actual cause.
    y = np.asarray(y)
    results = []
    all_subsets = [s for m in sizes for s in combinations(labels, m)]
    for n, subset in enumerate(all_subsets, 1):
        X = subset_matrix(subset, planet_mags, galaxy_mags)
        y_use, X_use = (y, X) if search_rows is None else (y[search_rows], X[search_rows])
        stats = score_subset(X_use, y_use, k_candidates=k_candidates,
                             reg_covar=reg_covar, n_seeds=n_seeds,
                             test_size=test_size, n_init=n_init)
        stats.update(subset=subset, n_filters=len(subset), n_colours=len(subset) - 1)
        results.append(stats)
        if progress and (n % progress == 0 or n == len(all_subsets)):
            print(f"  {n}/{len(all_subsets)} subsets")
    results.sort(key=lambda r: r["mean"], reverse=True)
    return results


def holdout_split(y, holdout_frac=0.3, random_state=0):
    """Split object indices once into (search pool, holdout), stratified.

    Called once, before any subset is looked at. Everything the search does
    happens inside the pool.
    """
    idx = np.arange(len(y))
    pool, holdout = train_test_split(
        idx, test_size=holdout_frac, stratify=y, random_state=random_state)
    return np.sort(pool), np.sort(holdout)


def evaluate_on_holdout(subset, planet_mags, galaxy_mags, y, pool_rows,
                        holdout_rows, *, k_candidates, reg_covar, n_init=10,
                        random_state=0):
    """Fit one subset on the whole search pool, score it once on the holdout.

    This is the only number that has not been selected on, and therefore the
    only one worth reporting.
    """
    X = subset_matrix(subset, planet_mags, galaxy_mags)
    k, clf, at_edge = select_k(X[pool_rows], y[pool_rows], k_candidates,
                               reg_covar, n_init=n_init, random_state=random_state)
    return {"subset": tuple(subset), "k": k, "k_at_edge": at_edge,
            "holdout_score": float(clf.score(X[holdout_rows], y[holdout_rows]))}


def winners_curse(n_candidates, sigma_per_split, n_seeds, trials=2000, seed=0):
    """Expected inflation of the best-of-N mean score, if all N were equal.

    A yardstick for reading the search table: differences at or below this are
    not distinguishable from picking the luckiest of N noisy draws.
    """
    rng = np.random.default_rng(seed)
    best = rng.normal(0.0, sigma_per_split,
                      size=(trials, n_candidates, n_seeds)).mean(axis=2).max(axis=1)
    return float(best.mean())
