"""Cached per-filter magnitudes for the planet and galaxy populations.

Why this exists
---------------
update_planet_flux_and_magnitude walks one planet and one filter at a time, and
galaxy_color_color_data_k15 walks one redshift and one filter at a time. Both
therefore cost O(n_filters) -- but each filter is computed independently of the
others, and colors are only differences of the results. So the expensive work
scales with the number of FILTERS, while the thing being studied is the number
of filter SUBSETS. Computing magnitudes once and slicing turns a subset sweep
from "recompute everything per subset" into pure arithmetic:

    table = planet_magnitudes("sonora-bobcat", ALL_FILTERS, **sample)
    # every subset below is now free
    for subset in combinations(ALL_FILTERS, 4):
        colors = add_color_columns(table, [filter_label(f) for f in subset])

Adding a filter later costs only that filter, not a full recompute.

The cache key
-------------
The dangerous failure here is silent: reuse magnitudes computed against a
DIFFERENT planet sample and nothing raises, the arrays are still the right
shape, and every downstream number is quietly wrong -- the same shape of bug as
the config seam that dwarfhunt.paths guards. So the key carries every input that
changes the sample (model tag, num_samples, radius_range, distance, rng seed,
deny-list), the full key is stored inside the cache file, and load_or_compute
re-checks it on read rather than trusting the filename hash. A mismatch raises.

Filters are deliberately NOT part of the key -- they are the incremental axis.
Two runs asking for different filter sets against the same sample share one
cache entry and each contributes the columns it computed.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
from species.read.read_model import ReadModel

from . import paths
from .galaxies import (check_filters_fit_k15_templates,
                       check_filters_fit_swire_templates, load_k15_data,
                       load_swire_data, redshift_data, redshift_k15_data,
                       synth_mags, translate_k15_L_v_to_f_lambda,
                       DEFAULT_REDSHIFTS)
from .planets import (filter_label, generate_planet_arrays,
                      update_planet_flux_and_magnitude)

CACHE_VERSION = 1


def default_cache_dir():
    """cache/photometry/, under the repo cache directory paths already owns."""
    return paths.cache_file("photometry")


def _key_digest(key):
    """Stable short hash of a key dict, for the filename only.

    Never the source of truth: the full key is written into the file and
    re-verified on load, so a hash collision surfaces as an error rather than as
    silently mismatched magnitudes.
    """
    blob = json.dumps(key, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:16]


def _normalise_seed(rng):
    """Reduce a seed to something the cache key can actually carry.

    `_key_digest` serialises with `default=str`, which quietly turns anything
    non-JSON into a string, and that breaks the key in two different ways:

    - `np.int64(2)` is stored as the string "2" but requested as np.int64(2), so
      `_load_entry`'s equality check fails and every read raises the alarming
      "written for a different sample" error even though the sample is fine.
    - a `Generator` stringifies with its memory address, so the filename hash
      differs every process. The cache never hits, a FRESH RANDOM SAMPLE is
      drawn on every run, and the numbers move with no signal at all -- exactly
      what the module docstring above promises cannot happen.

    So integers are normalised to one representation, and anything that cannot
    identify a sample is refused here rather than allowed to poison the key.
    """
    if rng is None:
        raise ValueError(
            "rng=None does not identify a sample, so the cache key cannot "
            "distinguish two different draws. Pass an integer seed.")

    if isinstance(rng, (int, np.integer)) and not isinstance(rng, bool):
        return int(rng)

    raise TypeError(
        f"rng must be an integer seed, got {type(rng).__name__}. Pass the seed "
        "itself, not a Generator: a Generator stringifies with its memory "
        "address, so every run would hash to a different cache entry and "
        "silently resample the population.")


def _load_entry(path, key):
    if not path.exists():
        return None
    with np.load(path, allow_pickle=False) as data:
        stored_key = json.loads(str(data["__key__"].item()))
        if stored_key != key:
            raise ValueError(
                f"cache file {path.name} was written for a different sample.\n"
                f"  stored:    {stored_key}\n"
                f"  requested: {key}\n"
                "Refusing to reuse it. Delete the file, or pass a cache_dir of "
                "your own.")
        return {k: data[k] for k in data.files if k != "__key__"}


def _save_entry(path, key, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(columns)
    payload["__key__"] = np.array(json.dumps(key, sort_keys=True, default=str))
    np.savez(path, **payload)


def planet_magnitudes(model_tag, filter_names, *, num_samples, radius_range,
                      rng, distance=10, deny=None, cache_dir=None,
                      refresh=False, verbose=True):
    """Planet sample plus abs_mag_/flux_ columns for every filter in `filter_names`.

    Returns the same flat dict shape update_planet_flux_and_magnitude produces,
    so it drops straight into add_color_columns / color_color_matrix.

    The sample itself (teff, logg, feh, radius, distance) is cached alongside the
    magnitudes and reused verbatim. Regenerating it from the seed would give the
    same draw today, but reusing the stored arrays means a filter added months
    later is guaranteed to describe the same objects even if the sampler's
    internals change.

    `rng` is a required integer seed, not optional and not a Generator -- see
    _normalise_seed. An unseeded draw cannot be named by the cache key, so
    `refresh=True` would swap the population out from under an unchanged key.
    """
    cache_dir = Path(cache_dir) if cache_dir else default_cache_dir()
    key = {
        "version": CACHE_VERSION,
        "kind": "planets",
        "model": model_tag,
        "num_samples": int(num_samples),
        "radius_range": [float(radius_range[0]), float(radius_range[1])],
        "distance": float(distance),
        "rng": _normalise_seed(rng),
        "deny": _key_digest(deny) if deny else None,
    }
    path = cache_dir / f"planets_{model_tag}_{_key_digest(key)}.npz"

    columns = None if refresh else _load_entry(path, key)
    model = None

    if columns is None:
        model = ReadModel(model_tag)
        sample = generate_planet_arrays(
            model, radius_range=radius_range, distance=distance,
            num_samples=num_samples, deny=deny, rng=rng)
        columns = {k: np.asarray(v) for k, v in sample.items()}

    wanted = [n for n in filter_names
              if f"abs_mag_{filter_label(n)}" not in columns]
    if wanted:
        if verbose:
            print(f"computing planet magnitudes for {len(wanted)} filter(s): "
                  f"{', '.join(filter_label(n) for n in wanted)}")
        model = model if model is not None else ReadModel(model_tag)
        # Feed the cached sample back in, so new filters describe the same
        # objects as the columns already stored.
        sample = {k: v for k, v in columns.items() if not k.startswith(("abs_mag_", "flux_"))}
        updated = update_planet_flux_and_magnitude(model, sample, wanted)
        for k, v in updated.items():
            columns[k] = np.asarray(v)
        _save_entry(path, key, columns)
    elif verbose:
        print(f"planet magnitudes: cache hit ({path.name})")

    return {k: v for k, v in columns.items()}


def galaxy_magnitudes(template, filter_names, *, redshifts=DEFAULT_REDSHIFTS,
                      cache_dir=None, refresh=False, verbose=True):
    """Redshift grid plus a mag_{label} column per filter for one K15 template.

    Colors are not included -- build them with planets.color_pairs so both
    populations go through the same primitive and the "A - B" keys line up.
    """
    cache_dir = Path(cache_dir) if cache_dir else default_cache_dir()
    template = Path(template)
    z = np.asarray(redshifts, dtype=float)
    key = {
        "version": CACHE_VERSION,
        "kind": "galaxies",
        "template": template.name,
        "redshifts": [float(z.min()), float(z.max()), int(z.size)],
    }
    path = cache_dir / f"galaxy_{template.stem}_{_key_digest(key)}.npz"

    # Coverage is a property of (filter, template, redshifts) and does not
    # depend on the cache, so check before either branch -- otherwise a cache
    # hit would skip the guard and hand back a column of NaN.
    check_filters_fit_k15_templates(template, filter_names, z)

    columns = None if refresh else _load_entry(path, key)
    if columns is None:
        columns = {"redshift": z}

    wanted = [n for n in filter_names if f"mag_{filter_label(n)}" not in columns]
    if wanted:
        if verbose:
            print(f"computing galaxy magnitudes for {template.name}, "
                  f"{len(wanted)} filter(s): "
                  f"{', '.join(filter_label(n) for n in wanted)}")
        wavelengths, fluxes, _errors = load_k15_data(template)
        labels = [filter_label(n) for n in wanted]
        cols = {label: [] for label in labels}
        for z_value in z:
            red_wl, red_flux = redshift_k15_data(wavelengths, fluxes, z_value)
            translated = translate_k15_L_v_to_f_lambda(red_wl, red_flux)
            for label, (mag, _err) in zip(labels, synth_mags(red_wl, translated, wanted)):
                cols[label].append(mag)
        for label, values in cols.items():
            columns[f"mag_{label}"] = np.asarray(values, dtype=float)
        _save_entry(path, key, columns)
    elif verbose:
        print(f"galaxy magnitudes: cache hit ({path.name})")

    return {k: v for k, v in columns.items()}


def galaxy_magnitudes_swire(template, filter_names, *, redshifts=DEFAULT_REDSHIFTS,
                            cache_dir=None, refresh=False, verbose=True):
    """Redshift grid plus a mag_{label} column per filter for one SWIRE template.

    The SWIRE counterpart to galaxy_magnitudes: identical cache mechanics, same
    return shape, so both feed planets.color_pairs the same way. Two differences
    from the K15 path:

    - SWIRE .sed flux is already F_lambda (erg cm^-2 s^-1 A^-1, normalised at
      5500 A per Polletta+2007), so there is no L_nu -> f_lambda translation
      step. The absolute normalisation is arbitrary and that is fine here --
      every number downstream is a colour, and a global scale cancels in a
      magnitude difference.
    - Coverage is checked with check_filters_fit_swire_templates. SWIRE spans
      ~0.1-6000 um rest-frame so it will not fire for the 2MASS/WISE set, but
      the guard still runs: a warm cache must not be a way to smuggle a
      truncated bandpass past the check.

    A distinct cache "kind" keeps SWIRE and K15 entries for a same-stemmed file
    from ever colliding, even though the extensions (.sed vs .txt) already
    differ.

    Colors are not included -- build them with planets.color_pairs so both
    populations go through the same primitive and the "A - B" keys line up.
    """
    cache_dir = Path(cache_dir) if cache_dir else default_cache_dir()
    template = Path(template)
    z = np.asarray(redshifts, dtype=float)
    key = {
        "version": CACHE_VERSION,
        "kind": "galaxies-swire",
        "template": template.name,
        "redshifts": [float(z.min()), float(z.max()), int(z.size)],
    }
    path = cache_dir / f"galaxy_{template.stem}_{_key_digest(key)}.npz"

    # Same reasoning as galaxy_magnitudes: coverage is a property of
    # (filter, template, redshifts) and independent of the cache, so check
    # before either branch or a cache hit would skip the guard.
    check_filters_fit_swire_templates(template, filter_names, z)

    columns = None if refresh else _load_entry(path, key)
    if columns is None:
        columns = {"redshift": z}

    wanted = [n for n in filter_names if f"mag_{filter_label(n)}" not in columns]
    if wanted:
        if verbose:
            print(f"computing galaxy magnitudes for {template.name}, "
                  f"{len(wanted)} filter(s): "
                  f"{', '.join(filter_label(n) for n in wanted)}")
        wavelengths, fluxes = load_swire_data(template)
        labels = [filter_label(n) for n in wanted]
        cols = {label: [] for label in labels}
        for z_value in z:
            red_wl, red_flux = redshift_data(wavelengths, fluxes, z_value)
            for label, (mag, _err) in zip(labels, synth_mags(red_wl, red_flux, wanted)):
                cols[label].append(mag)
        for label, values in cols.items():
            columns[f"mag_{label}"] = np.asarray(values, dtype=float)
        _save_entry(path, key, columns)
    elif verbose:
        print(f"galaxy magnitudes: cache hit ({path.name})")

    return {k: v for k, v in columns.items()}
