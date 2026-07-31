import json
from pathlib import Path


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

# Data files (missing_grid_points.json, species_database.hdf5) live in the
# notebook directory one level up from this helpers/ package, not wherever
# the importing notebook's cwd happens to be.
DATA_DIR = Path(__file__).resolve().parent.parent
DENY_PATH = DATA_DIR / "missing_grid_points.json"
DB_PATH = DATA_DIR / "species_database.hdf5"


def scan_missing_grid_points(tag, db_path=DB_PATH):
    """Find every grid cell of `tag` that holds an all-zero spectrum.

    Returns (params, axes, mask), where `mask` is a boolean array over the grid
    axes that is True wherever the spectrum is missing.
    """
    with h5py.File(db_path, "r") as hdf:
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


def load_deny_list(tags, db_path=DB_PATH, path=DENY_PATH, rebuild=False):
    """Load the cached deny-list, scanning the database the first time."""
    if path.exists() and not rebuild:
        return json.loads(path.read_text())

    deny = {}

    for tag in tags:
        params, axes, mask = scan_missing_grid_points(tag, db_path)

        combos = [
            [float(axes[param][i]) for param, i in zip(params, idx)]
            for idx in np.argwhere(mask)
        ]

        deny[tag] = {
            "params": params,
            "denied_axis_values": denied_axis_values(params, axes, mask),
            "combos": combos,
        }

        print(
            f"{tag}: {len(combos)} missing grid points"
            f" -> drop {deny[tag]['denied_axis_values']}"
        )

    path.write_text(json.dumps(deny, indent=1))

    return deny

deny = load_deny_list(['sonora-elfowl-t', 'sonora-elfowl-y'])


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


def generate_planet_arrays(model: ReadModel, radius_range, distance=10, num_samples=200, deny=None):
    model_bounds = model.get_bounds()
    model_points = model.get_points()

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

    random_values = {key: np.random.choice(vals, size=num_samples) for key, vals in pool.items()}
    #random_values = {key: np.random.uniform(low=val[0], high=val[1], size=num_samples) for key, val in model_bounds.items()}
    random_values['radius'] = np.random.uniform(low=radius_range[0], high=radius_range[1], size=num_samples)
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

def update_planet_flux_and_magnitude(model: ReadModel, planet_arrays):
    flux_model_1065 = ReadModel(model.model, wavel_range=(.61, 14.9), filter_name="JWST/MIRI.F1065C")
    flux_model_1140 = ReadModel(model.model, wavel_range=(.61, 14.9), filter_name="JWST/MIRI.F1140C")
    # unable to get this value to do model constraints
    #flux_model_1550 = ReadModel(model.model, wavel_range=(.61, 14.9), filter_name="JWST/MIRI.F1550C")

    flux_1065, flux_1140 = [], []
    abs_mag_1065, abs_mag_1140 = [], []

    with skip_nearest_spec_check():
        planet_count = 0
        # get_flux/get_magnitude take one planet at a time, so walk the rows
        for planet in iter_planets(planet_arrays):
            planet_flux_1065, planet_mag_1065 = flux_and_abs_mag(flux_model_1065, planet)
            planet_flux_1140, planet_mag_1140 = flux_and_abs_mag(flux_model_1140, planet)
            # unable to get this value to do model constraints
            #flux_1550.append(flux_model_1550.get_flux(planet)[0])

            flux_1065.append(planet_flux_1065)
            flux_1140.append(planet_flux_1140)
            abs_mag_1065.append(planet_mag_1065)
            abs_mag_1140.append(planet_mag_1140)
            planet_count += 1

    out = dict(planet_arrays)
    out['flux_1065C'] = np.array(flux_1065)
    out['flux_1140C'] = np.array(flux_1140)
    # unable to get this value to do model constraints
    #out['flux_1550C'] = np.array(flux_1550)
    out['abs_mag_1065'] = np.array(abs_mag_1065)
    out['abs_mag_1140'] = np.array(abs_mag_1140)

    out['F1065C - F1140C'] = out['abs_mag_1065'] - out['abs_mag_1140']

    dropped = int(np.count_nonzero(np.isnan(out['F1065C - F1140C'])))
    if dropped:
        print(f"{model.model}: {dropped} planets had no usable spectrum and are plotted as gaps")

    return out

def check_power_law(predicted_flux: np.ndarray, actual_flux: np.ndarray):
    # predicted_flux and actual_flux cover the same wavelengths, so this is how far
    # apart they are as a fraction of the real value -- the check of whether the
    # power law is a good fit. Broadcasts fine over the (n_planets, n_channels)
    # batches that get_planet_spectra/predict_power_law produce.
    residual = ((predicted_flux - actual_flux) / actual_flux)
    return residual


def get_planet_spectra(tag, planet_arrays, wavel_range=(0.61, 14.9), db_path=DB_PATH):
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
        axis_idx = [np.searchsorted(axes[param], planet_arrays[param]) for param in params]

        flux = np.empty((n_planets, hi - lo))
        for i in range(n_planets):
            flux[i] = group["flux"][
                axis_idx[0][i], axis_idx[1][i], axis_idx[2][i], axis_idx[3][i], axis_idx[4][i], lo:hi
            ]

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


def extrapolation_residuals(tag, planet_arrays, fit_range=(12.0, 14.0), test_range=(14.0, 14.9), db_path=DB_PATH):
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