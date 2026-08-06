"""Component-count-agnostic GMM classification, moved out of
michelson-gmm-exploration.ipynb.

The notebook's first pass fit a 2-component GaussianMixture and resolved the
arbitrary component naming with `1 - pred` -- a two-label flip that breaks the
moment n_components > 2 (component 3 maps to -2). That flip was standing in for
a more general idea: fit K unlabeled Gaussians, then ask which class each one
represents. This module does that generally, for any K and any number of
classes:

- fit a GaussianMixture on X only (labels never enter the fit)
- map each component to a class using training labels, via a class-balanced
  vote over soft responsibilities (predict_proba), not a hard-assignment
  majority vote
- expose predict / predict_proba / score / component_table on top

The class-balancing matters here specifically: dwarfs and galaxies are not
evenly represented (roughly 100 vs 270 in the current notebook), so a raw
vote would bias every mixed component toward "galaxy" regardless of what it
actually captures.

Nothing here assumes 2 feature columns or 2 classes -- component_class_map
and GMMClassifier work on whatever X and y are handed to them. The 2-D-only
assumption lives in color_color_plots.py's plotting helpers instead, since
that's where it's unavoidable.
"""

from typing import Literal

import numpy as np
from sklearn.mixture import GaussianMixture


def balanced_accuracy(y_true, y_pred, classes=None):
    """
    This function is used to compute the balanced accuracy score for the multiclass GMMs. 

    What this does, is pulls the correctly classified samples for each class, and the then divides it by the total number of samples in that class. 
    After it does that, it'll average those scores across all of the classes, giving you the balanced accuracy score. This works well for our case, 
    because we have classes that have more data points than others. This would still work if they were compeltely equal, but is generic. 

    Parameters
    ----------
    
    y_true: an array of the correct labels for the data points.
    y_pred: an array of the predicted labels for the data points. Created by the predict function.
    classes: an array of the unique classes in the data. If you don't pass this in, it'll unique the y_true array to get the actual classes you expect. 

    Returns
    -------

    balanced_accuracy: the balanced accuracy score for your GMM. 0.5 means that it's basically guessing, 1.0 means it's perfect. This value is a float. 

    """
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    if classes is None:
        classes = np.unique(y_true)

    recalls = []
    for c in classes:
        in_class = y_true == c
        n = in_class.sum()
        if n == 0:
            continue
        recalls.append(np.sum(in_class & (y_pred == c)) / n)
    return float(np.mean(recalls))


def component_class_map(resp, y_train, classes):
    """Assign each GMM component to the class it represents.

    Parameters
    ----------
    resp : ndarray, shape (n_train, n_components)
        Called the responsibility.
        This is the report card from the GMM training's data. Based on how much each point is assigned to each component,
        each row sums to 1.0 and each column sums to the total responsibility that component has for all of the training points.

    y_train : ndarray, shape (n_train,)
        This is the truth labels for the training data. This helps assign each component to a class.

    classes : ndarray
        Class values in a fixed order (index into this is the class index
        used everywhere else, e.g. mass's columns).

    Returns
    -------
    comp_to_class : ndarray, shape (n_components,)
        classes[j] that component k is mapped to.
        example: if comp_to_class[3] == 1, then component 3 is assigned to class 1.

    mass : ndarray, shape (n_components, n_classes)
        Class-balanced responsibility mass: for class c, the responsibility
        each component gets from c's training points, summed and divided by
        the count of the items in c. Balancing by class size keeps an under-represented class (e.g.
        the ~100 in class 1 against ~270 in class 2 here) from being outvoted in
        every mixed component just because it has fewer points.
    degenerate : ndarray of bool, shape (n_components,)
        True where a component drew ~no responsibility from any training
        point (mass row all zero). When this happens we shift it to the class with the most training points, 
        so it doesn't get left assigned to class 0
    """
    resp = np.asarray(resp)
    y_train = np.asarray(y_train)
    classes = np.asarray(classes)
    n_components = resp.shape[1]

    mass = np.zeros((n_components, len(classes)))
    for j, c in enumerate(classes):
        in_class = y_train == c
        n = in_class.sum()
        if n == 0:
            continue
        mass[:, j] = resp[in_class].sum(axis=0) / n

    # Which class has the most training points overall -- the fallback for
    # any component that gets no votes from anyone (see the loop below).
    class_counts = [np.sum(y_train == c) for c in classes]
    majority_class = classes[np.argmax(class_counts)]

    comp_to_class = np.empty(n_components, dtype=classes.dtype)
    degenerate = np.zeros(n_components, dtype=bool)

    for k in range(n_components):
        component_mass = mass[k]  # this component's vote share for each class
        got_any_votes = np.any(component_mass > 0)

        if got_any_votes:
            best_class_index = np.argmax(component_mass)
            comp_to_class[k] = classes[best_class_index]
        else:
            # No training point voted for this component at all -- argmax on
            # an all-zero row would silently pick class 0, so fall back to
            # the majority class instead and flag it as degenerate.
            comp_to_class[k] = majority_class
            degenerate[k] = True

    return comp_to_class, mass, degenerate


class GMMClassifier:
    """A GaussianMixture with an arbitrary number of components, wrapped so
    it behaves like a classifier.

    Fit is unsupervised, labels are used only afterward, to map
    components to classes. n_components can be larger than the number of
    classes -- the intended use here, since a single Gaussian per class
    can't trace a curved sequence (the dwarf T_eff track, the galaxy
    redshift tracks) and more components are what let the fit tile a curve
    with several small blobs instead of one badly-oriented ellipse.
    """

    def __init__(self, n_components=2, covariance_type: Literal['full', 'tied', 'diag', 'spherical'] = "full", n_init=10,
                 random_state=0):
        # n_init=10 (best-of-10-restarts by likelihood) is what actually
        # stabilizes the fit across seeds -- see Log.md 2026-08-05. Earlier
        # notes there attributed that stability to n_components; it was
        # n_init. Keeping the same default here so that behavior carries
        # forward unchanged for any K.
        self.n_components = n_components
        self.covariance_type : Literal['full', 'tied', 'diag', 'spherical'] = covariance_type
        self.n_init = n_init
        self.random_state = random_state

    def fit(self, X, y):
        X = np.asarray(X)
        y = np.asarray(y)

        self.gm_ = GaussianMixture(
            n_components=self.n_components,
            covariance_type=self.covariance_type,
            random_state=self.random_state,
            n_init=self.n_init,
        ).fit(X)

        self.classes_ = np.unique(y)
        resp = self.gm_.predict_proba(X)
        self.comp_to_class_, self.mass_, self.degenerate_ = component_class_map(
            resp, y, self.classes_)
        return self

    def predict(self, X):
        return self.comp_to_class_[self.gm_.predict(np.asarray(X))]

    def predict_proba(self, X):
        """Per-class probabilities, summing to 1 across class for each row of X. 

        Basically how likely each point is to be in each class. 

        Unlike predict (a hard argmax over components), this stays usable
        for anything that wants a soft score -- an ROC curve, a purity cut
        on the "which dwarfs fail" cell, etc.
        """
        resp = self.gm_.predict_proba(np.asarray(X))
        proba = np.zeros((resp.shape[0], len(self.classes_)))
        for j in range(len(self.classes_)):
            cols = np.where(self.comp_to_class_ == self.classes_[j])[0]
            if len(cols):
                proba[:, j] = resp[:, cols].sum(axis=1)
        return proba

    def score(self, X, y):
        return balanced_accuracy(y, self.predict(X), classes=self.classes_)

    def bic(self, X):
        return self.gm_.bic(np.asarray(X))

    def aic(self, X):
        return self.gm_.aic(np.asarray(X))

    def component_table(self, X_train=None, y_train=None):
        """One row per component: assigned class, class-balanced mass, raw
        counts, vote margin, mixture weight, and whether it was degenerate.

        Shows *which* components are ambiguous or mixed rather than only the
        aggregate score -- pass the same X_train/y_train used to fit() to
        get raw counts alongside the balanced mass.
        """
        sorted_mass = np.sort(self.mass_, axis=1)
        margin = sorted_mass[:, -1] - (sorted_mass[:, -2] if self.mass_.shape[1] > 1 else 0)

        raw_counts = None
        if X_train is not None and y_train is not None:
            hard = self.gm_.predict(np.asarray(X_train))
            y_train = np.asarray(y_train)
            raw_counts = np.zeros((self.n_components, len(self.classes_)), dtype=int)
            for k in range(self.n_components):
                in_comp = hard == k
                for j, c in enumerate(self.classes_):
                    raw_counts[k, j] = np.sum(in_comp & (y_train == c))

        rows = []
        for k in range(self.n_components):
            row = {
                "component": k,
                "assigned_class": self.comp_to_class_[k],
                "weight": self.gm_.weights_[k],
                "margin": margin[k],
                "degenerate": bool(self.degenerate_[k]),
            }
            for j, c in enumerate(self.classes_):
                row[f"mass_class_{c}"] = self.mass_[k, j]
                if raw_counts is not None:
                    row[f"n_class_{c}"] = raw_counts[k, j]
            rows.append(row)
        return rows


def fit_and_score(X_train, y_train, X_test, y_test, n_components=2, seed=0,
                   **kwargs):
    """Fit a GMMClassifier and score it on held-out data.

    Kept as a module-level function with the same (score, pred) return shape
    the notebook's seed-stability loop already depends on -- now a thin
    wrapper over GMMClassifier, with n_components exposed (default 2, so
    existing 2-component calls are unchanged).
    """
    clf = GMMClassifier(n_components=n_components, random_state=seed, **kwargs)
    clf.fit(X_train, y_train)
    pred = clf.predict(X_test)
    return balanced_accuracy(y_test, pred, classes=clf.classes_), pred
