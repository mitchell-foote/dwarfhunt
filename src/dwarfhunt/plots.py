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

import numpy as np
from matplotlib import colormaps
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


def _covariances_2d(gm):
    """Normalize any GaussianMixture covariance_type to (n_components, 2, 2).

    sklearn stores covariances in four different shapes depending on
    covariance_type ("full": (K,2,2), "tied": (2,2), "diag": (K,2), "spherical":
    (K,)); the ellipse math below needs a full 2x2 per component regardless of
    which one was used to fit, otherwise it only happens to work for "full".
    """
    cov = gm.covariances_
    k = gm.n_components
    if gm.covariance_type == "full":
        return cov
    if gm.covariance_type == "tied":
        return np.broadcast_to(cov, (k, *cov.shape)).copy()
    if gm.covariance_type == "diag":
        out = np.zeros((k, 2, 2))
        out[:, 0, 0] = cov[:, 0]
        out[:, 1, 1] = cov[:, 1]
        return out
    if gm.covariance_type == "spherical":
        out = np.zeros((k, 2, 2))
        out[:, 0, 0] = cov
        out[:, 1, 1] = cov
        return out
    raise ValueError(f"unrecognized covariance_type: {gm.covariance_type!r}")


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
                            n_std=2, linewidth=1.4):
    """Draw each GMM component as a covariance ellipse.

    Built for the K-components-per-curve GMMClassifier in gmm_classify.py --
    this is how you actually see whether the fitted components are tiling a
    population's curve or cutting across it, rather than just reading off an
    aggregate accuracy number.

    Parameters
    ----------
    ax : matplotlib Axes
    gm : sklearn.mixture.GaussianMixture (already fit, 2 features)
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

    Returns
    -------
    list of Ellipse patches added to ax, one per component.
    """
    covariances = _covariances_2d(gm)
    means = gm.means_
    weights = gm.weights_
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
                          point_alpha=0.7):
    """Shade the color-color plane by predicted class and scatter the real
    points on top.

    clf needs a .predict(X) -> class array method (GMMClassifier from
    gmm_classify.py fits this directly). Only meaningful in 2 features, since
    it shades the axes' own plane.

    Builds the meshgrid from the *current* axis limits via get_xlim/get_ylim
    rather than assuming ascending order, so this still works after
    style_panel has inverted the y axis for the magnitude convention.

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
    if X.shape[1] != 2:
        raise ValueError(f"plot_decision_regions needs 2 feature columns, got {X.shape[1]}")

    classes = np.unique(y)
    if class_colors is None:
        class_colors = class_color_map(classes)

    x_lo, x_hi = sorted(ax.get_xlim())
    y_lo, y_hi = sorted(ax.get_ylim())
    xx, yy = np.meshgrid(np.linspace(x_lo, x_hi, resolution),
                         np.linspace(y_lo, y_hi, resolution))
    grid_pred = clf.predict(np.column_stack((xx.ravel(), yy.ravel())))

    class_index = {c: i for i, c in enumerate(classes)}
    zz = np.vectorize(class_index.get)(grid_pred).reshape(xx.shape)

    region_cmap = ListedColormap([class_colors[c] for c in classes])
    ax.contourf(xx, yy, zz, levels=np.arange(len(classes) + 1) - 0.5,
               cmap=region_cmap, alpha=0.25, zorder=0)
    for c in classes:
        mask = y == c
        ax.scatter(X[mask, 0], X[mask, 1], color=class_colors[c], s=16,
                  alpha=point_alpha, edgecolors="0.25", linewidths=0.3,
                  zorder=3, label=str(c))
    return ax
