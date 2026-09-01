import json


from contextlib import contextmanager
from itertools import product
import species.read.read_model as read_model_module
from species.read.read_model import ReadModel

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from species.core import constants
from species.read.read_filter import ReadFilter

from . import paths

# --- deny-list of grid points with no spectrum --------------------------------
# species stores an all-zero spectrum for any grid point whose model file is
# absent and that linear interpolation could not recover. That is the
# "Could not interpolate N grid points so storing zeros instead" warning that
# add_model prints. Sampling one of those points gives a band flux of exactly
# 0.0, and get_magnitude then dies inside math.log10(0.0) with
# "math domain error". get_flux does not raise on them, it just hands back 0.0,
# so they quietly poison the flux columns as well.
#
# For Elf Owl every missing point sits at logg = 3.0: the distribution ships
# only 12 logg_3.00 spectra (all feh=1.0, co=0.458, logkzz=8.0) out of 1440,
# and logg=3.0 is on the edge of the axis so griddata could not extrapolate
# the rest. Nothing at logg >= 3.25 is missing.

# Data locations come from dwarfhunt.paths, which reads the same
# species_config.ini that species itself reads. The previous
# `Path(__file__).parent.parent` form silently repointed the moment this module
# moved, and gave no error when the file it named was absent -- it just rebuilt
# the cache against nothing.
#
# These are functions rather than module constants because the database path is
# only knowable after dwarfhunt.init() has chosen a config, and importing this
# module must not depend on that having happened yet.


def default_db_path():
    """The configured database. Resolved per call, never cached at import."""
    return paths.database_path()


def default_deny_path():
    """The deny-list cache under the repo cache directory."""
    return paths.cache_file("missing_grid_points.json")


def scan_missing_grid_points(tag, db_path=None):
    """Find every grid cell of `tag` that holds an all-zero spectrum.

    Returns (params, axes, mask), where `mask` is a boolean array over the grid
    axes that is True wherever the spectrum is missing.
    """
    db_path = db_path if db_path is not None else default_db_path()

    with h5py.File(db_path, "r") as hdf:
        if f"models/{tag}" not in hdf:
            # h5py's own KeyError names the missing tag but not the alternatives,
            # and a typo'd tag is the usual cause.
            available = sorted(hdf["models"]) if "models" in hdf else []
            raise KeyError(
                f"{db_path} has no model {tag!r}. Available: {available}")

        group = hdf[f"models/{tag}"]
        params = [group.attrs[f"parameter{i}"] for i in range(group.attrs["n_param"])]
        axes = {param: group[param][:] for param in params}

        flux = group["flux"]
        mask = np.zeros(flux.shape[: len(params)], dtype=bool)

        # One teff slice at a time (~68 MB) instead of the whole ~818 MB array.
        for i in range(flux.shape[0]):
            mask[i] = ~np.any(flux[i] != 0.0, axis=-1)

    return params, axes, mask


def denied_axis_values(params, axes, mask):
    """Smallest set of axis values to drop so no missing cell stays reachable.

    Greedy set cover: repeatedly drop whichever (parameter, value) accounts for
    the most still-missing cells. Because the sampler draws each axis
    independently, removing these values from the pool makes every possible
    draw valid by construction.
    """
    remaining = np.argwhere(mask)
    dropped = {}

    while len(remaining):
        best_axis, best_idx, best_count = None, None, -1

        for axis, param in enumerate(params):
            for idx in np.unique(remaining[:, axis]):
                count = int(np.count_nonzero(remaining[:, axis] == idx))
                if count > best_count:
                    best_axis, best_idx, best_count = axis, idx, count

        best_param = params[best_axis]
        dropped.setdefault(best_param, []).append(float(axes[best_param][best_idx]))
        remaining = remaining[remaining[:, best_axis] != best_idx]

    return dropped


def _scan_one(tag, db_path):
    """Scan a single tag and return its deny-list entry."""
    params, axes, mask = scan_missing_grid_points(tag, db_path)

    combos = [
        [float(axes[param][i]) for param, i in zip(params, idx)]
        for idx in np.argwhere(mask)
    ]

    entry = {
        "params": params,
        "denied_axis_values": denied_axis_values(params, axes, mask),
        "combos": combos,
    }

    print(
        f"{tag}: {len(combos)} missing grid points"
        f" -> drop {entry['denied_axis_values']}"
    )

    return entry


def load_deny_list(tags, db_path=None, path=None, rebuild=False):
    """Deny-list entries for `tags`, scanning the database for anything missing.

    The cache file is a {tag: entry} map and is checked TAG BY TAG, never
    wholesale. An earlier version returned `json.loads(path.read_text())` as soon
    as the file existed, whatever had been asked for -- so once the file held the
    two Elf Owl tags, asking for any other tag handed back the Elf Owl entries
    and no rescan happened. generate_planet_arrays then found no entry for its
    model, applied no denial at all, and could sample a grid point with no
    spectrum; the only symptom was NaN magnitudes appearing much later.

    So: tags already in the file are reused, tags that are not are scanned and
    merged in, and the return value holds exactly `tags` -- a missing key is
    then impossible rather than silent. Other tags stay in the file, since
    scanning one costs a pass over the whole flux array.
    """
    db_path = db_path if db_path is not None else default_db_path()
    path = path if path is not None else default_deny_path()

    cached = json.loads(path.read_text()) if path.exists() else {}

    missing = [t for t in tags if rebuild or t not in cached]
    if missing:
        for tag in missing:
            cached[tag] = _scan_one(tag, db_path)
        path.write_text(json.dumps(cached, indent=1))

    return {tag: cached[tag] for tag in tags}


# The Elf Owl grids are the ones with missing spectra, so they are what the
# cached deny-list covers by default.
DEFAULT_DENY_TAGS = ("sonora-elfowl-t", "sonora-elfowl-y")

_DENY_CACHE = {}


def deny_list(tags=DEFAULT_DENY_TAGS, rebuild=False):
    """The deny-list, built on first use and memoised for the session.

    This used to run at module scope. That made a bare `import` of this module
    open the multi-GB database and write a JSON cache as a side effect -- and it
    fired transitively, so importing the galaxy module paid for it too. Worse,
    it ran before any caller could choose a config, so it read whichever
    database the cwd implied.
    """
    key = tuple(tags)

    if rebuild or key not in _DENY_CACHE:
        _DENY_CACHE[key] = load_deny_list(list(tags), rebuild=rebuild)

    return _DENY_CACHE[key]


@contextmanager
def skip_nearest_spec_check():
    """Turn off species' per-call missing-spectrum check inside this block.

    This is what made each magnitude take about a minute. ReadModel.get_model
    calls check_nearest_spec on every single call, and that helper builds a
    fresh ReadModel and calls get_data once per corner of the 5-D grid cell, so
    2**5 = 32 times. get_data does

        flux = np.array(hdf5_file[f"models/{self.model}/flux"])
        flux = flux[..., self.wl_index]

    i.e. it reads the whole ~818 MB flux dataset off disk and only then keeps
    the filter's handful of channels. So one get_flux or get_magnitude costs
    about 26 GB of I/O, and one planet costs four of those.

    Measured on this grid: 10.7 s per call with the check, 0.0026 s without,
    and the returned flux and magnitude are bit-identical. All the check does
    is warn that a nearest grid point has no spectrum, which the deny-list
    above already guarantees cannot happen for the points we sample, so there
    is nothing for it to find. That warning is also the only thing it emits, so
    this silences the noise as a side effect.
    """
    original = read_model_module.check_nearest_spec
    read_model_module.check_nearest_spec = lambda *args, **kwargs: None
    try:
        yield
    finally:
        read_model_module.check_nearest_spec = original


def generate_planet_arrays(model: ReadModel, radius_range, distance=10, num_samples=200, deny=None, rng=None):
    model_bounds = model.get_bounds()
    model_points = model.get_points()
    generator = np.random.default_rng(rng) if rng is not None else np.random.default_rng()
    # Strip the values with no spectrum out of the pool before sampling, so
    # every draw is valid by construction. This throws away the handful of real
    # spectra that happen to share a denied value (Elf Owl ships 12 logg=3.00
    # files), but reaching one needs three other parameters to coincide as
    # well, so independent per-axis sampling would essentially never hit them.
    model_deny = (deny or {}).get(model.model, {})
    denied = model_deny.get('denied_axis_values', {})
    pool = {}
    for key in model_bounds.keys():
        values = model_points[key]
        keep = values[~np.isin(values, denied.get(key, []))]
        if len(keep) < len(values):
            print(f"{model.model}: dropping {key} = {sorted(set(values) - set(keep))} (no spectra)")
        pool[key] = keep

    # nothing denied should still be reachable from the filtered pool
    denied_combos = {tuple(combo) for combo in model_deny.get('combos', [])}
    if denied_combos:
        params = model_deny['params']
        reachable = product(*(pool[param] for param in params))
        assert not any(tuple(float(v) for v in combo) in denied_combos for combo in reachable), \
            "deny-list does not cover the whole sampling pool"

    random_values = {key: generator.choice(vals, size=num_samples) for key, vals in pool.items()}
    random_values['radius'] = generator.uniform(low=radius_range[0], high=radius_range[1], size=num_samples)
    random_values['distance'] = np.full(num_samples, distance)
    return random_values

def iter_planets(planet_arrays):
      n = len(next(iter(planet_arrays.values())))
      for i in range(n):
          yield {key: float(val[i]) for key, val in planet_arrays.items()}

def flux_and_abs_mag(reader: ReadModel, planet):
    # get_flux returns (flux, uncertainty) and the uncertainty is always None
    band_flux = reader.get_flux(planet)[0]

    if not np.isfinite(band_flux) or band_flux <= 0.0: # type: ignore
        # a missing spectrum stored as zeros, get_magnitude would raise
        # "math domain error" out of log10, so leave a gap instead
        return np.nan, np.nan

    # get_magnitude returns (apparent, absolute)
    return band_flux, reader.get_magnitude(planet)[1] # type: ignore


def filter_label(filter_name):
    """"JWST/MIRI.F1065C" -> "F1065C". Shared by every flux/mag/color column name."""
    return filter_name.rsplit(".", 1)[-1]


def filter_mean_wavelength(filter_name):
    """Mean wavelength (um) of a filter's bandpass, from its transmission curve.

    The sort key for put_filters_in_wavelength_order. Read from species rather
    than parsed out of the filter name: the name is only a rough guide and is
    actively misleading for broad bands. "WISE/WISE.W3" is nominally the 12 um
    band and sorts after "F1140C" alphabetically *and* numerically, but its
    bandpass spans 7.2 - 18.4 um with a mean of 12.8 um, which places it between
    F1140C (11.3) and F1550C (15.5) -- not at the red end where the raw name
    suggests.
    """
    return float(ReadFilter(filter_name).mean_wavelength())


def put_filters_in_wavelength_order(filter_names):
    """Sort filters blue to red by mean wavelength.

    color_pairs documents its convention as "bluer minus redder, short
    wavelength first", but it can only honour that if the labels it is handed
    are already in wavelength order -- it pairs by position, not by physics. Run
    the filter list through here before deriving labels and colors and the
    convention holds automatically, whatever order the filters were typed in.

    This matters more as filters are added: with an unsorted list, two different
    filter subsets can produce adjacent-colour bases ordered differently, which
    makes their results awkward to compare for no real reason.

    Returns a tuple, so it can be used as a module-level constant without the
    aliasing risk a list would carry.
    """
    return tuple(sorted(filter_names, key=filter_mean_wavelength))


def assert_wavelength_ordered(filter_names):
    """Raise unless `filter_names` runs blue to red. Returns them unchanged.

    The companion check to put_filters_in_wavelength_order. Everything
    downstream that builds adjacent-pair colours -- color_pairs' "bluer minus
    redder" convention, search.colour_names, and the n-1 colour basis the
    subset sweep is built on -- pairs filters BY POSITION and cannot tell
    whether the positions mean anything. Hand it an unordered list and the
    colours are still computed, still named, and still fitted; some are simply
    the negative of what their name says.

    That is a sign flip with no exception, which is why this is a hard check
    rather than a warning. It takes full species filter names because that is
    the only form a mean wavelength can be read from -- a bare label like
    "F1065C" carries no wavelength, and the name is not a reliable guide anyway
    (see filter_mean_wavelength for the WISE.W3 case).
    """
    names = list(filter_names)
    means = [filter_mean_wavelength(n) for n in names]

    out_of_order = [i for i in range(len(means) - 1) if means[i] > means[i + 1]]
    if out_of_order:
        shown = ", ".join(f"{n} ({m:.2f} um)" for n, m in zip(names, means))
        raise ValueError(
            f"filters are not in wavelength order: {shown}. Adjacent-pair "
            "colours are built by position, so this would silently flip the "
            "sign of the colours across "
            f"{', '.join(f'{names[i]}/{names[i + 1]}' for i in out_of_order)}. "
            "Run the list through put_filters_in_wavelength_order first."
        )

    return names


def model_wavel_range(tag, db_path=None):
    """Wavelength coverage (um) actually stored for `tag` in the database.

    Reads just the wavelength dataset (a few thousand floats), not the ~818 MB flux
    array, so this is cheap enough to call as a precondition check. This is the ground
    truth for what a grid can support -- ReadModel's own wavel_range argument is not:
    it's silently overwritten by the filter's own range whenever filter_name is set
    (see check_filters_fit_model), so it can't be used to detect a mismatch itself.
    """
    db_path = db_path if db_path is not None else default_db_path()

    with h5py.File(db_path, "r") as hdf:
        wl = hdf[f"models/{tag}/wavelength"]
        return float(wl[0]), float(wl[-1])


def check_filters_fit_model(tag, filter_names, db_path=None):
    """Raise ValueError for any filter whose bandpass falls outside `tag`'s grid.

    Elf Owl (0.61-14.9um) can't support JWST/MIRI.F1550C (14.912-16.219um) for exactly
    this reason -- the grid stops 12nm short of where the filter starts. Bobcat
    (0.61-17.0um) covers it fine. Failing fast here beats the "math domain error" that
    would otherwise surface deep inside get_magnitude.
    """
    model_lo, model_hi = model_wavel_range(tag, db_path=db_path)

    for name in filter_names:
        filt_lo, filt_hi = ReadFilter(name).wavelength_range()
        if filt_lo < model_lo or filt_hi > model_hi:
            raise ValueError(
                f"{tag} covers {model_lo:.3f}-{model_hi:.3f} um but {name} needs "
                f"{filt_lo:.3f}-{filt_hi:.3f} um. Re-add the model with a wider "
                "wavel_range, or drop this filter for this model."
            )


def update_planet_flux_and_magnitude(
    model: ReadModel,
    planet_arrays,
    filter_names=("JWST/MIRI.F1065C", "JWST/MIRI.F1140C"),
):
    """Add per-filter flux_{label} and abs_mag_{label} columns to `planet_arrays`.

    filter_names can be any length -- pass all three MIRI coronagraph filters for
    Bobcat, or just the two that fit Elf Owl's narrower grid (see
    check_filters_fit_model). Color columns are not computed here; call
    add_color_columns on the result.
    """
    check_filters_fit_model(model.model, filter_names)

    # wavel_range is not passed here: ReadModel overwrites it from the filter's own
    # range whenever filter_name is set, so passing one would be dead code.
    readers = [ReadModel(model.model, filter_name=name) for name in filter_names]
    labels = [filter_label(name) for name in filter_names]

    flux_cols = {label: [] for label in labels}
    mag_cols = {label: [] for label in labels}

    with skip_nearest_spec_check():
        # get_flux/get_magnitude take one planet at a time, so walk the rows
        for planet in iter_planets(planet_arrays):
            for reader, label in zip(readers, labels):
                planet_flux, planet_mag = flux_and_abs_mag(reader, planet)
                flux_cols[label].append(planet_flux)
                mag_cols[label].append(planet_mag)

    out = dict(planet_arrays)
    for label in labels:
        out[f"flux_{label}"] = np.array(flux_cols[label])
        out[f"abs_mag_{label}"] = np.array(mag_cols[label])

    dropped = int(np.count_nonzero([np.isnan(out[f"abs_mag_{label}"]) for label in labels]))
    if dropped:
        print(f"{model.model}: {dropped} planet/filter magnitudes had no usable spectrum and are gaps")

    return out


def color_pairs(mags_by_filter, order=None):
    """Every pairwise magnitude difference, bluer minus redder, short wavelength first.

    mags_by_filter : dict, filter label -> magnitude (scalar or ndarray)
        e.g. {"F1065C": m1, "F1140C": m2, "F1550C": m3}. Works unchanged on scalars
        (one galaxy template at one redshift) or ndarrays (a batch of planets),
        since it's pure subtraction.
    order : sequence of labels, optional
        Pins the pair ordering. Defaults to mags_by_filter's insertion order, i.e.
        the order filter_names was given in.

    Returns
    -------
    dict, "{blue} - {red}" -> difference, one entry per i < j pair.
    """
    labels = list(order) if order is not None else list(mags_by_filter.keys())
    return {
        f"{labels[i]} - {labels[j]}": mags_by_filter[labels[i]] - mags_by_filter[labels[j]]
        for i in range(len(labels))
        for j in range(i + 1, len(labels))
    }


def add_color_columns(planet_data, filter_names=None):
    """Return a copy of `planet_data` with every pairwise color column added.

    planet_data : dict
        Output of update_planet_flux_and_magnitude (or anything with abs_mag_{label}
        keys). teff/logg/radius/flux_* etc. all pass through unchanged, so this stays
        one flat table per model -- what plotting and the eventual galaxy merge want.
    filter_names : list of filter labels, optional
        Restrict/order which abs_mag_{label} columns are combined, e.g.
        ["F1065C", "F1140C", "F1550C"]. Defaults to every abs_mag_* column present,
        in insertion order.

    Colors are differences, so this is safe to compare against the galaxy notebook's
    apparent-magnitude colors even though planets here use absolute magnitude -- the
    distance modulus cancels either way.
    """
    prefix = "abs_mag_"
    if filter_names is None:
        filter_names = [key[len(prefix):] for key in planet_data if key.startswith(prefix)]

    mags_by_filter = {label: planet_data[f"{prefix}{label}"] for label in filter_names}

    out = dict(planet_data)
    out.update(color_pairs(mags_by_filter, order=filter_names))
    return out


def color_color_matrix(planet_data, color_names=None):
    """Stack this model's color columns into one (n_planets, n_colors) array.

    add_color_columns leaves colors as separate dict entries ("F1065C - F1140C",
    "F1140C - F1550C", ...), convenient for name-based lookup but awkward for
    scatter code that wants to index by column position -- see
    michelson-galaxy-graph.ipynb's plot_color_color, which does exactly that
    against a per-galaxy [F1065-F1140, F1140-F1550, F1065-F1550] tuple. This is the
    batched, planet-array equivalent: one row per planet instead of one tuple per
    galaxy/redshift, indexed the same way.

    Parameters
    ----------
    planet_data : dict
        Output of add_color_columns -- needs one array per color, keyed "A - B".
    color_names : list of str, optional
        Which color columns to stack, and in what order, e.g.
        ["F1065C - F1140C", "F1140C - F1550C", "F1065C - F1550C"] to match the
        galaxy notebook's axis order. Defaults to every "A - B" key present, in
        insertion order.

    Returns
    -------
    color_names : list[str]
        Resolved column order -- color_names[j] labels colors[:, j].
    colors : ndarray, shape (n_planets, n_colors)

    Example
    -------
    >>> names, colors = color_color_matrix(planet_colors)
    >>> x = colors[:, names.index("F1065C - F1140C")]
    >>> y = colors[:, names.index("F1065C - F1550C")]
    >>> ax.scatter(x, y, c=planet_colors["teff"], cmap="viridis")
    """
    if color_names is None:
        color_names = [key for key in planet_data if " - " in key]

    colors = np.column_stack([np.atleast_1d(planet_data[name]) for name in color_names])
    return color_names, colors


def check_power_law(predicted_flux: np.ndarray, actual_flux: np.ndarray):
    # predicted_flux and actual_flux cover the same wavelengths, so this is how far
    # apart they are as a fraction of the real value -- the check of whether the
    # power law is a good fit. Broadcasts fine over the (n_planets, n_channels)
    # batches that get_planet_spectra/predict_power_law produce.
    residual = ((predicted_flux - actual_flux) / actual_flux)
    return residual


def _axis_indices(axis, values, tag, param, rtol=1e-9):
    """Locate each value on a grid axis, refusing anything not actually on it.

    np.searchsorted alone is not enough. It returns an *insertion point*, not a
    match, so a value that misses a node returns the index of a neighbouring one
    and the caller silently gets the wrong grid point's spectrum -- no exception,
    no warning, just a subtly wrong number. It can also return len(axis) for a
    value past the top end, which then raises a confusing out-of-bounds error far
    from the cause.

    This function is only valid because generate_planet_arrays draws every
    parameter with np.random.choice over the grid's own axis values, so every
    requested value is expected to sit exactly on a node. If that assumption ever
    breaks, this raises instead of quietly lying.
    """
    axis = np.asarray(axis)
    values = np.asarray(values)

    idx = np.searchsorted(axis, values)
    idx = np.clip(idx, 0, len(axis) - 1)

    # searchsorted can land one past the true node for a value equal to it,
    # depending on side; check the neighbour too before declaring a miss.
    left = np.clip(idx - 1, 0, len(axis) - 1)
    use_left = np.abs(axis[left] - values) < np.abs(axis[idx] - values)
    idx = np.where(use_left, left, idx)

    off_node = ~np.isclose(axis[idx], values, rtol=rtol, atol=0.0)
    if np.any(off_node):
        bad = np.unique(values[off_node])[:5]
        raise ValueError(
            f"{tag}: {int(off_node.sum())} value(s) for {param!r} are not on the "
            f"model grid, e.g. {bad.tolist()}. get_planet_spectra indexes the grid "
            "directly and cannot interpolate; sample with generate_planet_arrays, "
            "or go through species.ReadModel.get_model instead."
        )

    return idx


def get_planet_spectra(tag, planet_arrays, wavel_range=(0.61, 14.9), db_path=None):
    """Read a batch of planet spectra straight out of the HDF5 grid, no species.

    generate_planet_arrays samples every parameter with np.random.choice against
    the grid's own axis values, so each planet lands exactly on a grid node --
    there's no interpolation to do. That makes it safe to skip
    species.ReadModel.get_data, which reads the *entire* ~818 MB flux dataset off
    disk on every single call regardless of how many planets you ask for (see
    skip_nearest_spec_check's docstring for where that number comes from). Indexing
    the array directly instead takes ~0.04s for 200 planets instead of tens of
    minutes.

    Parameters
    ----------
    tag : str
        Model tag, e.g. "sonora-elfowl-t".
    planet_arrays : dict
        Output of generate_planet_arrays: one array per grid parameter (all grid
        values, since these are sampled from the grid) plus 'radius' and 'distance'.
    wavel_range : (float, float)
        Wavelength window (um) to return, inclusive.
    db_path : Path

    Returns
    -------
    wl : ndarray, shape (n_channels,)
        Wavelengths (um) within `wavel_range`, shared by every planet.
    flux : ndarray, shape (n_planets, n_channels)
        Flux density at 10 pc scaled per planet by (radius/distance)**2, using the
        same constants and formula as ReadModel.get_data so the two are
        bit-identical at a shared grid point.
    """
    db_path = db_path if db_path is not None else default_db_path()

    with h5py.File(db_path, "r") as hdf:
        group = hdf[f"models/{tag}"]
        params = [group.attrs[f"parameter{i}"] for i in range(group.attrs["n_param"])]
        axes = {param: group[param][:] for param in params}

        wl_full = group["wavelength"][:]
        channel_idx = np.where((wl_full >= wavel_range[0]) & (wl_full <= wavel_range[1]))[0]
        lo, hi = channel_idx[0], channel_idx[-1] + 1
        wl = wl_full[lo:hi]

        n_planets = len(next(iter(planet_arrays.values())))
        # each parameter's grid value -> its index along that axis, per planet
        axis_idx = [_axis_indices(axes[param], planet_arrays[param], tag, param)
                    for param in params]

        flux = np.empty((n_planets, hi - lo))
        for i in range(n_planets):
            # Build the index tuple from however many parameters this grid has.
            # Hardcoding five positions only ever worked for Elf Owl: bobcat is
            # 3-parameter (flux is 4-D) and diamondback is 4-parameter, and both
            # raised "list index out of range" here.
            flux[i] = group["flux"][tuple(ax[i] for ax in axis_idx) + (slice(lo, hi),)]

    # (radius/distance)**2 scaling -- same formula and constants ReadModel.get_data
    # uses (read_model.py, "Apply (radius/distance)^2 scaling"), so a direct
    # HDF5 read and a species get_data call agree bit-for-bit on the same planet.
    scaling = (planet_arrays["radius"] * constants.R_JUP) ** 2 / (
        planet_arrays["distance"] * constants.PARSEC
    ) ** 2
    flux = flux * scaling[:, None]

    return wl, flux


def get_planet_spectra_via_species(tag, planet_arrays, wavel_range=(0.61, 14.9)):
    """Same interface and output as get_planet_spectra, but through species.ReadModel.get_data.

    Exists purely to check the fast direct-HDF5 path agrees with species -- not for
    routine use. get_data reads the *entire* ~818 MB flux dataset off disk on every
    call regardless of how many planets you ask for, so this costs ~18s for 200
    planets on one grid, vs get_planet_spectra's ~0.03s. It does not need
    skip_nearest_spec_check: that guards check_nearest_spec, which only runs inside
    get_model, not get_data.

    Returns
    -------
    wl : ndarray, shape (n_channels,)
    flux : ndarray, shape (n_planets, n_channels)
    """
    model = ReadModel(tag, wavel_range=wavel_range)
    n_planets = len(next(iter(planet_arrays.values())))

    # species pads a few extra points past wavel_range on each side (for
    # filter-profile resampling); trim back to the plain inclusive range
    # get_planet_spectra uses so the two line up channel-for-channel.
    wl_padded, _ = model.wavelength_points()
    keep = (wl_padded >= wavel_range[0]) & (wl_padded <= wavel_range[1])
    wl = wl_padded[keep]

    flux = np.empty((n_planets, wl.size))
    for i, planet in enumerate(iter_planets(planet_arrays)):
        flux[i] = model.get_data(planet).flux[keep]

    return wl, flux


def fit_power_law(wl, flux, fit_range):
    """Fit flux = 10**intercept * wl**slope per planet, in log-log space.

    Vectorized across planets: np.polyfit accepts a 2-D `y`, one column per
    planet, and returns coefficients shaped (2, n_planets).

    Parameters
    ----------
    wl : ndarray, shape (n_channels,)
    flux : ndarray, shape (n_planets, n_channels)
    fit_range : (float, float)
        Wavelength window (um) to fit on.

    Returns
    -------
    slope, intercept : ndarray, shape (n_planets,)
        NaN for any planet whose fit window contains a non-positive or
        non-finite flux value (log10 would be undefined for it).
    """
    mask = (wl >= fit_range[0]) & (wl <= fit_range[1])
    fit_flux = flux[:, mask]

    n_planets = flux.shape[0]
    slope = np.full(n_planets, np.nan)
    intercept = np.full(n_planets, np.nan)

    valid = np.all(np.isfinite(fit_flux) & (fit_flux > 0), axis=1)
    if np.any(valid):
        coef = np.polyfit(np.log10(wl[mask]), np.log10(fit_flux[valid]).T, 1)
        slope[valid], intercept[valid] = coef[0], coef[1]

    return slope, intercept


def predict_power_law(wl_target, slope, intercept):
    """Evaluate a fitted power law at `wl_target` for every planet.

    Returns shape (n_planets, len(wl_target)); NaN rows propagate from a planet
    that fit_power_law couldn't fit.
    """
    return 10 ** (intercept[:, None] + slope[:, None] * np.log10(wl_target)[None, :])


def extrapolation_residuals(tag, planet_arrays, fit_range=(12.0, 14.0), test_range=(14.0, 14.9), db_path=None):
    """Fit a power law on `fit_range` and check it against real flux on `test_range`.

    This is the validation step: unlike the 14.9-16.6um region species can't cover
    at all, `test_range` sits inside the grid, so `actual` here is real model flux
    rather than another extrapolation -- it's the one place we can check the method
    against ground truth before trusting it past the grid edge.

    Returns
    -------
    dict with:
        tag, wl_test : the test-range wavelengths
        residual : ndarray (n_planets, n_test), see check_power_law
        slope, intercept : ndarray (n_planets,), the fitted power law per planet
        teff : ndarray (n_planets,), for coloring plots by temperature
    """
    db_path = db_path if db_path is not None else default_db_path()

    wl, flux = get_planet_spectra(tag, planet_arrays, wavel_range=(fit_range[0], test_range[1]), db_path=db_path)

    slope, intercept = fit_power_law(wl, flux, fit_range)

    test_mask = (wl >= test_range[0]) & (wl <= test_range[1])
    wl_test = wl[test_mask]
    actual = flux[:, test_mask]
    predicted = predict_power_law(wl_test, slope, intercept)

    return {
        "tag": tag,
        "wl_test": wl_test,
        "residual": check_power_law(predicted, actual),
        "slope": slope,
        "intercept": intercept,
        "teff": np.asarray(planet_arrays["teff"], dtype=float),
    }


def residual_by_teff(result):
    """Bin a result dict's residual by teff instead of by wavelength.

    Works on any result dict shaped like extrapolation_residuals' output (needs just
    "residual" and "teff"), so the same call works for the Bobcat proxy check too.
    generate_planet_arrays samples teff from the grid's own axis values, so grouping
    by the exact value recovers the grid's actual temperature steps rather than
    needing arbitrary bin edges.

    For each teff, all (planet, wavelength) residuals at planets with that teff are
    pooled together -- same flattening extrapolation_residuals' aggregate stats use,
    just narrowed to one temperature at a time. Percentiles are of the signed
    residual (so median always falls between p16 and p84); p84 |residual| is kept
    separate as "how bad does a typical bad channel get."

    Returns a list of dicts, one per distinct teff, sorted low to high.
    """
    teff = result["teff"]
    residual_pct = result["residual"] * 100

    rows = []
    for t in np.unique(teff):
        r = residual_pct[teff == t]
        rows.append({
            "teff": float(t),
            "n_planets": int(np.count_nonzero(teff == t)),
            "p16 residual (%)": float(np.nanpercentile(r, 16)),
            "median residual (%)": float(np.nanmedian(r)),
            "p84 residual (%)": float(np.nanpercentile(r, 84)),
            "p84 |residual| (%)": float(np.nanpercentile(np.abs(r), 84)),
        })

    return rows


def wavelength_flux_power_law(wl, flux, fit_range, wl_ext_stop, n_ext=200):
    """Extend `flux` past its native coverage with a fitted power law.

    This is the production counterpart to extrapolation_residuals: same
    fit_power_law/predict_power_law core, but extending past the grid's edge
    (e.g. to 16.6um for F1550C) instead of into a region we can check against
    real data. Only trust this as far as extrapolation_residuals showed the
    method holds up.

    Parameters
    ----------
    wl : ndarray, shape (n_channels,)
    flux : ndarray, shape (n_planets, n_channels)
    fit_range : (float, float)
        Wavelength window (um) to fit the power law on.
    wl_ext_stop : float
        Wavelength (um) to extend out to.
    n_ext : int
        Number of extension points.

    Returns
    -------
    wl_full, flux_full : the input concatenated with the extrapolated tail.
    """
    slope, intercept = fit_power_law(wl, flux, fit_range)

    wl_ext = np.linspace(wl[-1], wl_ext_stop, n_ext)
    flux_ext = predict_power_law(wl_ext, slope, intercept)

    wl_full = np.concatenate([wl, wl_ext])
    flux_full = np.concatenate([flux, flux_ext], axis=1)
    return wl_full, flux_full


def plot_extrapolation_residuals(result, title=None, ax=None):
    """Plot one line per planet of extrapolation residual vs wavelength.

    Colored by teff (viridis: perceptually uniform and colorblind-safe) so a
    temperature-correlated bias -- the kind already documented for the
    Rayleigh-Jeans extension -- would show up as banding rather than noise.
    Median and 16th-84th percentile band summarize the ~200 individual planet
    lines, which are too dense to read on their own.
    """
    wl = result["wl_test"]
    residual_pct = result["residual"] * 100
    teff = result["teff"]

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 5))

    norm = Normalize(vmin=teff.min(), vmax=teff.max())
    cmap = mpl.colormaps["viridis"] # type: ignore

    for row, t in zip(residual_pct, teff):
        ax.plot(wl, row, color=cmap(norm(t)), alpha=0.15, linewidth=0.6)

    median = np.nanmedian(residual_pct, axis=0)
    p16 = np.nanpercentile(residual_pct, 16, axis=0)
    p84 = np.nanpercentile(residual_pct, 84, axis=0)
    ax.fill_between(wl, p16, p84, color="0.3", alpha=0.2, label="16th-84th pct")
    ax.plot(wl, median, color="0.1", linewidth=1.8, label="median")
    ax.axhline(0, color="0.5", linewidth=1, linestyle="--")

    # (predicted - actual) / actual is formally unbounded wherever a deep absorption
    # trough puts `actual` near zero -- a real feature of the spectrum, not a bug --
    # so a handful of channels can blow out to 4-5 digit percentages and, left to
    # autoscale, crush every other channel to a flat line. Clip the view to a robust
    # range instead; the individual planet lines still reach past it, so the extreme
    # channels are visible as lines running off the top/bottom, just not dictating the
    # scale everyone else is squeezed into.
    lo, hi = np.nanpercentile(residual_pct, [1, 99])
    pad = 0.15 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("Teff (K)")

    ax.set_xlabel("Wavelength (µm)")
    ax.set_ylabel("Residual (%)  [(predicted − actual) / actual]")
    ax.set_title(title or f"{result['tag']}: power-law extrapolation residual")
    ax.legend(loc="upper left", frameon=False)

    if fig is not None:
        fig.tight_layout()
        return fig, ax
    return ax