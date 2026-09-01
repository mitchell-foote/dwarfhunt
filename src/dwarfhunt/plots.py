"""Plotting for the brown-dwarf vs. galaxy color-color figures.

Split out of generate_galaxy_data.py so that module stays about *producing*
magnitudes and colors, and this one is about drawing them. Nothing here
touches species, templates, or filters -- it only consumes the dicts and
tuple-lists the data helpers hand back, so it imports no other helper module.

Styled to match the Michelson figure:
- y axes run downward, the usual magnitude convention (bright/blue at top)
- galaxy tracks are continuous lines colour-graded by redshift, red -> yellow
- brown dwarfs are a dark -> light blue ramp in effective temperature
- two horizontal colourbars sit above the panels
- tracks are labelled directly rather than through a shape legend
"""

from itertools import combinations

import numpy as np
from matplotlib import colormaps
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse


def _clip_cmap(name, lo, hi, label):
    """Take a slice out of a built-in ramp.

    Both ends of the two ramps used here run to near-white, which disappears
    against the page. Clipping keeps the full ramp readable on white.
    """
    return LinearSegmentedColormap.from_list(
        label, colormaps[name](np.linspace(lo, hi, 256)))


# Dark navy (cool) -> light blue (hot). Clipped short of white so the hottest
# dwarfs still register as marks; the dark marker edge does the rest.
DWARF_CMAP = _clip_cmap("Blues_r", 0.0, 0.82, "dwarf_blues")
# Red (low z) -> yellow (high z), as in the reference figure.
GALAXY_CMAP = colormaps["autumn"]


def track_xy(data, x_idx, y_idx, x_label, y_label):
    """Pull (x, y, redshift) out of either shape the data helpers produce.

    - dict from galaxy_color_color_data* -- keyed by "A - B" label, so
      x_idx/y_idx are ignored.
    - list of (z, mag_array, color_color_data) tuples from
      get_full_redshift_mag_loop* -- x_idx/y_idx index into each row's
      color_color_data (0=F1065C-F1140C, 1=F1140C-F1550C, 2=F1065C-F1550C).
    """
    if isinstance(data, dict):
        return (np.asarray(data[x_label]), np.asarray(data[y_label]),
                np.asarray(data["redshift"]))
    z_values = np.asarray([item[0] for item in data])
    color_color = [item[2] for item in data]
    return (np.asarray([row[x_idx] for row in color_color]),
            np.asarray([row[y_idx] for row in color_color]),
            z_values)


def plot_color_color(data, title, ax, color, x_idx, y_idx, x_label, y_label):
    """Scatter one template's color-color track in a single flat color.

    The original plotting path, kept for michelson-galaxy-graph.ipynb's style
    of figure. plot_galaxy_track below is the redshift-encoded replacement.

    Feeding the dict shape into the x_idx/y_idx path used to fail with
    "string index out of range": iterating a dict walks its string keys, not
    rows, so item[2] grabbed a single character instead of a color value.
    """
    x, y, _z_values = track_xy(data, x_idx, y_idx, x_label, y_label)
    ax.plot(x, y, color=color, alpha=0.3, label=title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend()


def plot_dwarfs(ax, x, y, teff, norm, cmap=DWARF_CMAP):
    """Scatter the brown dwarfs, coloured by effective temperature.

    The thin dark edge is what keeps the hot (pale) end of the ramp visible
    against white. Returns the handle a colorbar needs.
    """
    return ax.scatter(x, y, c=teff, cmap=cmap, norm=norm, s=22,
                      edgecolors="0.25", linewidths=0.4, zorder=3)


def plot_galaxy_track(data, ax, x_idx, y_idx, x_label, y_label, norm,
                      cmap=GALAXY_CMAP, label=None, linewidth=1.6,
                      label_offset=(6, 6)):
    """Draw one galaxy template's redshift track as a colour-graded line.

    Built as a LineCollection so colour varies *along* the line rather than
    per marker -- the track is a continuous path through colour-colour space,
    and drawing it as discrete points hid that. Each segment takes the mean
    redshift of its two endpoints.

    A dot marks the z-minimum end so you can tell which way the track runs,
    and `label` is written there directly. Direct labels beat a legend here:
    the three tracks share one colormap, so a legend swatch could not tell
    them apart anyway.

    Returns the handle a colorbar needs.
    """
    x, y, z_values = track_xy(data, x_idx, y_idx, x_label, y_label)

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    line = LineCollection(segments, cmap=cmap, norm=norm, zorder=2)
    line.set_array(0.5 * (z_values[:-1] + z_values[1:]))
    line.set_linewidth(linewidth)
    ax.add_collection(line)

    start = int(np.argmin(z_values))
    ax.scatter(x[start], y[start], c=[z_values[start]], cmap=cmap, norm=norm,
               s=18, zorder=4)
    if label:
        ax.annotate(label, (x[start], y[start]), textcoords="offset points",
                    xytext=label_offset, fontsize=8.5, color="0.25")

    # LineCollection does not drive autoscaling on its own.
    ax.autoscale_view()
    return line


def style_panel(ax, x_label, y_label):
    """Axis labels plus the downward y axis the magnitude convention wants.

    Call after everything is drawn: invert_yaxis flips the *current* limits,
    so inverting before the data is in place would be undone by autoscaling.
    """
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.margins(0.06)
    ax.autoscale_view()
    if not ax.yaxis_inverted():
        ax.invert_yaxis()


def add_figure_scales(fig, axes, dwarf_handle, galaxy_handle,
                      legend_ax=-1, top=0.82, gap=0.05):
    """Two horizontal colourbars above the panels, plus a small legend.

    Placed with explicit figure coordinates rather than by handing matplotlib
    the axes list -- that route sizes each bar by how many panels it is
    attached to, which would make the two scales different widths.

    `gap` is the space between the top of the panels and the bottom of the
    bars; the tick labels and title stack *above* each bar, so the bars need
    headroom to the figure edge, not below.
    """
    fig.subplots_adjust(top=top)
    cax_galaxy = fig.add_axes([0.11, top + gap, 0.34, 0.026])
    cax_dwarf = fig.add_axes([0.56, top + gap, 0.34, 0.026])

    for handle, cax, title in ((galaxy_handle, cax_galaxy, "Galaxy Redshift"),
                               (dwarf_handle, cax_dwarf, "Brown Dwarf $T_{eff}$")):
        cbar = fig.colorbar(handle, cax=cax, orientation="horizontal")
        cax.xaxis.set_ticks_position("top")
        cax.xaxis.set_label_position("top")
        cbar.set_label(title, labelpad=8)
        cbar.outline.set_visible(False)

    handles = [
        Line2D([], [], color=GALAXY_CMAP(0.0), linewidth=1.6, label="Galaxy track"),
        Line2D([], [], marker="o", linestyle="none", color=DWARF_CMAP(0.35),
               markeredgecolor="0.25", markeredgewidth=0.4, markersize=5,
               label="Brown dwarf"),
    ]
    axes[legend_ax].legend(handles=handles, loc="upper right", frameon=True,
                           framealpha=0.9, fontsize=9)
    return cax_galaxy, cax_dwarf


def _covariances_full(gm):
    """Normalize any GaussianMixture covariance_type to (n_components, d, d).

    sklearn stores covariances in four different shapes depending on
    covariance_type ("full": (K,d,d), "tied": (d,d), "diag": (K,d), "spherical":
    (K,)); the ellipse and marginalization math below needs a dense d x d per
    component regardless of which one was used to fit, otherwise it only happens
    to work for "full".
    """
    cov = np.asarray(gm.covariances_)
    k = gm.n_components
    d = np.asarray(gm.means_).shape[1]
    diag = np.arange(d)

    if gm.covariance_type == "full":
        return cov
    if gm.covariance_type == "tied":
        return np.broadcast_to(cov, (k, d, d)).copy()
    if gm.covariance_type == "diag":
        out = np.zeros((k, d, d))
        out[:, diag, diag] = cov
        return out
    if gm.covariance_type == "spherical":
        out = np.zeros((k, d, d))
        out[:, diag, diag] = cov[:, None]
        return out
    raise ValueError(f"unrecognized covariance_type: {gm.covariance_type!r}")


def marginal_mixture(gm, dims):
    """The fitted mixture restricted to `dims`, as (weights, means, covariances).

    Marginalizing a Gaussian mixture onto a subset of its coordinates is exact,
    not an approximation: integrating the other coordinates out leaves a mixture
    with the *same* weights, the sub-vector of each mean, and the corresponding
    sub-block of each covariance. That is what makes a 2-color panel of a
    d-color fit an honest picture rather than a cartoon -- it is precisely the
    model the classifier would have if it had only ever seen those two colors.

    What it is *not* is a slice through the d-dimensional decision boundary; see
    plot_decision_regions's `mode` for that distinction.
    """
    dims = list(dims)
    weights = np.asarray(gm.weights_)
    means = np.asarray(gm.means_)[:, dims]
    covariances = _covariances_full(gm)[:, dims][:, :, dims]
    return weights, means, covariances


def _weighted_log_prob(points, weights, means, covariances):
    """log(w_k * N(x; mu_k, Sigma_k)), shape (n_points, n_components).

    The same quantity sklearn's GaussianMixture argmaxes in predict(), rebuilt
    here because the marginal mixture from marginal_mixture() is a bare
    (weights, means, covariances) triple with no fitted estimator to ask.
    """
    points = np.asarray(points, dtype=float)
    n, d = points.shape
    out = np.empty((n, len(weights)))

    with np.errstate(divide="ignore"):
        log_weights = np.log(weights)

    for k in range(len(weights)):
        cov = covariances[k]
        sign, logdet = np.linalg.slogdet(cov)
        if sign <= 0:
            # A component can go singular in a marginal projection even when the
            # d-dimensional fit is healthy. Nudge the diagonal rather than
            # emitting a nan that would silently blank part of the panel.
            cov = cov + np.eye(d) * 1e-12
            sign, logdet = np.linalg.slogdet(cov)
        delta = points - means[k]
        maha = np.einsum("ij,jk,ik->i", delta, np.linalg.inv(cov), delta)
        out[:, k] = log_weights[k] - 0.5 * (maha + logdet + d * np.log(2 * np.pi))
    return out


def class_color_map(classes):
    """Fixed class -> color assignment, shared by every class-colored plot
    here (component ellipses, decision regions).

    Both call this with the same `classes` array when the caller doesn't pin
    its own class_colors, so ellipses, the decision-region background, and
    the scattered points always agree on which color means which class --
    passing a colormap name straight to contourf does not give that
    guarantee (see plot_decision_regions).
    """
    cycle = colormaps["tab10"].colors
    return {c: cycle[i % len(cycle)] for i, c in enumerate(sorted(classes))}


def plot_component_ellipses(ax, gm, comp_to_class=None, class_colors=None,
                            n_std=2, linewidth=1.4, dims=(0, 1)):
    """Draw each GMM component as a covariance ellipse.

    Built for the K-components-per-curve GMMClassifier in gmm_classify.py --
    this is how you actually see whether the fitted components are tiling a
    population's curve or cutting across it, rather than just reading off an
    aggregate accuracy number.

    Parameters
    ----------
    ax : matplotlib Axes
    gm : sklearn.mixture.GaussianMixture (already fit, any number of features)
    comp_to_class : ndarray, shape (n_components,), optional
        Component -> class assignment from gmm_classify.component_class_map
        or GMMClassifier.comp_to_class_. Colors ellipses by class instead of
        by component index when given.
    class_colors : dict, class -> color, optional
        Only used when comp_to_class is given. Defaults to matplotlib's
        default color cycle, keyed by sorted class value.
    n_std : float
        Ellipse radius in standard deviations (2 ~ 95% contour for a 2-D
        Gaussian).
    dims : (int, int)
        Which two feature columns this panel is drawing. For a fit with more
        than 2 features the ellipse is the exact marginal of the component onto
        those two columns (see marginal_mixture), i.e. the shadow the ellipsoid
        casts on this plane -- so it is drawn at full extent, not narrowed the
        way a slice through the ellipsoid would be.

    Returns
    -------
    list of Ellipse patches added to ax, one per component.
    """
    weights, means, covariances = marginal_mixture(gm, dims)
    max_weight = weights.max() if len(weights) else 1.0

    if comp_to_class is not None and class_colors is None:
        class_colors = class_color_map(set(comp_to_class))

    patches = []
    for k in range(gm.n_components):
        vals, vecs = np.linalg.eigh(covariances[k])
        vals = np.clip(vals, 0, None)  # guard tiny negative round-off
        order = vals.argsort()[::-1]
        vals, vecs = vals[order], vecs[:, order]
        angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
        width, height = 2 * n_std * np.sqrt(vals)

        color = "0.35" if comp_to_class is None else class_colors[comp_to_class[k]]
        ellipse = Ellipse(means[k], width, height, angle=angle,
                          facecolor=color, edgecolor=color,
                          alpha=0.15 + 0.35 * (weights[k] / max_weight),
                          linewidth=linewidth, zorder=2.5)
        ax.add_patch(ellipse)
        patches.append(ellipse)
    return patches


def plot_decision_regions(ax, clf, X, y, resolution=300, class_colors=None,
                          point_alpha=0.7, dims=(0, 1), mode="marginal",
                          slice_at=None):
    """Shade a color-color plane by predicted class and scatter the real
    points on top.

    clf needs a .predict(X) -> class array method (GMMClassifier from gmm.py
    fits this directly). `X` is always the *full* feature matrix; `dims` picks
    the two columns this panel draws.

    Builds the meshgrid from the *current* axis limits via get_xlim/get_ylim
    rather than assuming ascending order, so this still works after
    style_panel has inverted the y axis for the magnitude convention.

    Parameters
    ----------
    dims : (int, int)
        Feature columns for the x and y axis.
    mode : {"marginal", "slice"}
        Only consulted when X has more than 2 columns, because with exactly 2
        the panel *is* the feature space and the regions are the classifier's
        real ones. Beyond that, a plane cannot show a d-dimensional boundary
        and you have to say which 2-D question you are asking:

        - "marginal" (default): integrate the other colors out under the
          fitted mixture, then decide. Answers "what would this classifier do
          knowing only these two colors?". Uses every training point, so the
          shading is stable, but it is a projection -- two regions that look
          overlapping here may be cleanly separated by a color not on display.
        - "slice": hold the other colors fixed at `slice_at` and evaluate the
          real d-dimensional classifier there. Answers "where is the true
          boundary on this particular cut?". Genuinely a cross-section of the
          boundary, but it only describes that one cut, and most plotted points
          do not lie on it.

        Neither is the whole boundary. That is a property of d > 2, not a
        shortcoming of either choice.
    slice_at : array-like, shape (n_features,), optional
        Where to hold the off-panel columns for mode="slice". Defaults to the
        per-column median of X, i.e. the middle of the observed data.
    class_colors : dict, class -> color, optional
        Defaults to class_color_map(classes). Pass the same dict used for
        plot_component_ellipses's comp_to_class coloring to keep background,
        ellipses, and scattered points in agreement -- handing a named
        colormap straight to contourf does *not* give that for free: contourf
        normalizes across the colormap's full range, so with e.g. "tab10" (10
        entries) and only 2 classes it can pick two arbitrary, non-adjacent
        colors instead of the map's first two.
    """
    X = np.asarray(X)
    y = np.asarray(y)
    x_dim, y_dim = dims
    n_features = X.shape[1]

    classes = np.unique(y)
    if class_colors is None:
        class_colors = class_color_map(classes)

    # Scatter first, then read the limits. This function used to build the
    # meshgrid before anything had been drawn, so on a fresh axes get_xlim()
    # returned matplotlib's default (0, 1) and the shading came out as a
    # rectangle from 0 to 1 that ignored where the data actually sat. Plotting
    # the points first lets autoscale establish real limits; contourf goes
    # underneath afterward via zorder=0, so the draw order is unchanged.
    for c in classes:
        mask = y == c
        ax.scatter(X[mask, x_dim], X[mask, y_dim], color=class_colors[c], s=16,
                  alpha=point_alpha, edgecolors="0.25", linewidths=0.3,
                  zorder=3, label=str(c))

    x_lo, x_hi = sorted(ax.get_xlim())
    y_lo, y_hi = sorted(ax.get_ylim())
    # Pad beyond the current view so the shading still reaches the corners after
    # style_panel applies its margins, rather than leaving an unshaded border.
    pad_x = 0.12 * (x_hi - x_lo)
    pad_y = 0.12 * (y_hi - y_lo)
    xx, yy = np.meshgrid(np.linspace(x_lo - pad_x, x_hi + pad_x, resolution),
                         np.linspace(y_lo - pad_y, y_hi + pad_y, resolution))
    grid = np.column_stack((xx.ravel(), yy.ravel()))

    if n_features == 2:
        # The panel is the whole feature space: ask the classifier directly.
        grid_pred = clf.predict(grid)
    elif mode == "slice":
        if slice_at is None:
            slice_at = np.median(X, axis=0)
        full = np.tile(np.asarray(slice_at, dtype=float), (len(grid), 1))
        full[:, [x_dim, y_dim]] = grid
        grid_pred = clf.predict(full)
    elif mode == "marginal":
        if not hasattr(clf, "gm_"):
            raise ValueError(
                'mode="marginal" needs a GMMClassifier (it marginalizes clf.gm_); '
                'pass mode="slice" for a classifier that only exposes .predict')
        # Hard argmax over components, then map to class -- deliberately the
        # same rule GMMClassifier.predict uses, so the shading agrees with the
        # confusion matrix rather than quietly using a softer criterion.
        weights, means, covariances = marginal_mixture(clf.gm_, dims)
        log_prob = _weighted_log_prob(grid, weights, means, covariances)
        grid_pred = clf.comp_to_class_[log_prob.argmax(axis=1)]
    else:
        raise ValueError(f'mode must be "marginal" or "slice", got {mode!r}')

    # The shaded regions come from the CLASSIFIER, which knows every class it was
    # trained on -- not just the ones present in the `y` handed to this panel. So
    # the index has to span both. Building it from np.unique(y) alone made
    # np.vectorize(class_index.get) return None for any class the classifier
    # predicted but the plotted subset did not contain, and the only symptom was
    # "TypeError: int() argument must be ... not 'NoneType'" from inside
    # np.vectorize, which names nothing useful. Plotting one class's failures, or
    # any row subset, is enough to hit it.
    region_classes = list(classes)
    region_classes += [c for c in np.unique(grid_pred) if c not in region_classes]

    # Colour any extra class without disturbing what the caller pinned, so a
    # shared class_colors dict still means the same thing in every panel.
    if len(region_classes) > len(classes):
        fallback = class_color_map(region_classes)
        class_colors = {**{c: fallback[c] for c in region_classes}, **class_colors}

    class_index = {c: i for i, c in enumerate(region_classes)}
    zz = np.vectorize(class_index.get)(grid_pred).reshape(xx.shape)

    region_cmap = ListedColormap([class_colors[c] for c in region_classes])
    keep_x, keep_y = ax.get_xlim(), ax.get_ylim()
    ax.contourf(xx, yy, zz, levels=np.arange(len(region_classes) + 1) - 0.5,
               cmap=region_cmap, alpha=0.25, zorder=0)
    ax.set_xlim(keep_x)
    ax.set_ylim(keep_y)
    return ax


def plot_failures(ax, X, y_true, y_pred, cls, dims=(0, 1), c=None,
                  cmap=None, norm=None, marker="o", size=28,
                  highlight_color="red", subset=None):
    """Scatter one class's points on a 2-color projection, ringing the misses.

    The generic form of the "which dwarfs fail" / "which galaxies fail" panels:
    those differed only in which class they focused on and what they colored the
    points by (T_eff for dwarfs, redshift for galaxies), so both are this call
    with different `cls` and `c`.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features)
        Full feature matrix; `dims` picks the two columns drawn.
    y_true, y_pred : ndarray, shape (n_samples,)
    cls : the class value to draw. Rows of another class are ignored entirely.
    c : ndarray, shape (n_samples,), optional
        Per-point value for the color scale, in the *same row order as X* --
        it is masked to `cls` here, so pass the full-length array, not a
        pre-filtered one.
    subset : ndarray of bool, shape (n_samples,), optional
        Further restrict which rows are drawn, on top of the class mask. Lets
        one class be split across several calls that differ only in marker --
        e.g. one per K15 AGN-fraction template -- while redshift keeps the
        color channel to itself. The returned counts describe the drawn subset,
        not the whole class.

    Returns
    -------
    (handle, n_wrong, n_total) -- handle is what a colorbar needs.
    """
    X = np.asarray(X)
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    x_dim, y_dim = dims

    sel = y_true == cls
    if subset is not None:
        sel = sel & np.asarray(subset)
    x = X[sel, x_dim]
    y_vals = X[sel, y_dim]
    wrong = y_pred[sel] != cls

    scatter_kw = dict(marker=marker, s=size, edgecolors="0.25",
                      linewidths=0.4, zorder=3)
    if c is None:
        handle = ax.scatter(x, y_vals, **scatter_kw)
    else:
        handle = ax.scatter(x, y_vals, c=np.asarray(c)[sel], cmap=cmap,
                            norm=norm, **scatter_kw)

    ax.scatter(x[wrong], y_vals[wrong], facecolors="none",
               edgecolors=highlight_color, marker=marker, s=size * 3.4,
               linewidths=1.4, zorder=4)
    return handle, int(wrong.sum()), int(sel.sum())


def feature_pairs(n_features):
    """Every (i, j) column pair, i < j -- one panel per pair."""
    return list(combinations(range(n_features), 2))


def plot_decision_grid(clf, X, y, names, pairs=None, mode="marginal",
                       slice_at=None, class_colors=None, class_labels=None,
                       ncols=3, panel_size=4.4, resolution=200,
                       ellipses=True, suptitle=None):
    """One decision-region panel per pair of color axes.

    The replacement for a single hardcoded (names[0], names[1]) figure. With d
    colors there are C(d, 2) planes to look at and no principled reason to
    privilege one, so draw them all: 2 colors gives 1 panel (identical to the
    old figure), 3 gives 3, 4 gives 6.

    Every panel shares one `class_colors` mapping so a color means the same
    class throughout, and every panel is titled with its own axis pair -- with
    six lookalike panels, an unlabeled one is worse than no panel.

    Returns
    -------
    (fig, axes, pairs) -- axes[i] draws pairs[i], so a caller can add per-panel
    annotation afterward without recomputing the pairing.
    """
    X = np.asarray(X)
    y = np.asarray(y)
    if pairs is None:
        pairs = feature_pairs(X.shape[1])

    if class_colors is None:
        class_colors = class_color_map(np.unique(y))

    ncols = min(ncols, len(pairs))
    nrows = int(np.ceil(len(pairs) / ncols))
    # layout="constrained" rather than a later tight_layout(): the panels carry
    # long y labels ("F1140C - F1550C") and constrained layout measures them
    # before placing axes, instead of fitting them into a fixed grid afterward.
    # tight_layout also cannot account for a colorbar spanning several axes.
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(panel_size * ncols, panel_size * nrows),
                             squeeze=False, layout="constrained")
    axes = axes.ravel()

    for ax, (x_dim, y_dim) in zip(axes, pairs):
        plot_decision_regions(ax, clf, X, y, dims=(x_dim, y_dim), mode=mode,
                              slice_at=slice_at, class_colors=class_colors,
                              resolution=resolution)
        if ellipses and hasattr(clf, "gm_"):
            plot_component_ellipses(ax, clf.gm_, comp_to_class=clf.comp_to_class_,
                                    class_colors=class_colors,
                                    dims=(x_dim, y_dim))
        style_panel(ax, names[x_dim], names[y_dim])

    for ax in axes[len(pairs):]:
        ax.set_visible(False)

    if class_labels:
        handles = [Line2D([], [], marker="o", linestyle="none",
                          markerfacecolor=class_colors[c], markeredgecolor="0.25",
                          markersize=7, label=class_labels.get(c, str(c)))
                   for c in sorted(class_colors)]
        axes[0].legend(handles=handles, title="true class", fontsize=8.5)

    if suptitle:
        fig.suptitle(suptitle)
    return fig, axes[:len(pairs)], pairs
