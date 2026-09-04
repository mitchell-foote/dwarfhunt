"""Galaxy template helpers (SWIRE and Kirkpatrick+2015), moved out of
michelson-galaxy-graph.ipynb.

Mirrors planets.py's shape: load a template, redshift it, synthesize
per-filter magnitudes, then hand the magnitudes to that module's color_pairs so
planets and galaxies compute colors the same way.

Data only -- the drawing lives in color_color_plots.py.
"""

from itertools import combinations
from pathlib import Path

import numpy as np
from species.phot.syn_phot import SyntheticPhotometry
from species.read.read_filter import ReadFilter

from .planets import color_pairs, filter_label

DEFAULT_FILTERS = ("JWST/MIRI.F1065C", "JWST/MIRI.F1140C", "JWST/MIRI.F1550C")
DEFAULT_REDSHIFTS = np.linspace(0.5, 2, 90)


def load_swire_data(file):
    """
    Load the SWIRE data from a text file, based on https://www.iasf-milano.inaf.it/~polletta/templates/swire_templates.html

    The templates should be in the format of a text file with the following columns:
    1. Wavelength (in Angstroms)
    2. Flux (normalized at 5500 Angstroms in erg / cm^2 / s / Angstrom)

    Returns ( wavelengths in microns, fluxes in erg / cm^2 / s / Angstrom )
    """
    data = np.loadtxt(file, comments='#', usecols=(0, 1))
    wavelengths = data[:, 0] * 1e-4
    fluxes = data[:, 1]
    return wavelengths, fluxes


def redshift_data(wavelengths, fluxes, z):
    """
    Redshift the data by a given redshift z.

    Parameters:
    - wavelengths: array of wavelengths in microns
    - fluxes: array of fluxes in erg / cm^2 / s / Angstrom
    - z: redshift value

    Returns:
    - redshifted_wavelengths: array of redshifted wavelengths in microns
    - redshifted_fluxes: array of redshifted fluxes in erg / cm^2 / s / Angstrom
    """
    redshifted_wavelengths = wavelengths * (1 + z)
    redshifted_fluxes = fluxes
    return redshifted_wavelengths, redshifted_fluxes


def synth_mags(wavelength, flux, filter_names):
    """
    Calculate synthetic magnitudes for a given list of filters.

    spectrum_to_magnitude returns ((apparent, error), (absolute, error)); this keeps
    only the apparent-magnitude tuple per filter, same as the original notebook cell
    -- callers that just want the scalar apparent magnitude take item[0] of each
    entry (see get_full_redshift_mag_loop and galaxy_color_color_data below).
    """
    mag_list = []
    for filter_name in filter_names:
        filters_data = SyntheticPhotometry(filter_name).spectrum_to_magnitude(wavelength, flux)
        mag_list.append(filters_data[0])
    return mag_list


def get_color_color_data(mag_array):
    """
    Calculate color-color data from an array of magnitudes.

    Parameters:
    - mag_array: array of magnitudes in the order [F1065C, F1140C, F1550C]

    Returns:
    - color_color_data: array of color-color data in the order
      [F1065C-F1140C, F1140C-F1550C, F1065C-F1550C]

    Kept as-is (fixed positional order) because plot_color_color and the existing
    notebook cell index into it by position. galaxy_color_color_data below is the
    named, color_pairs-based equivalent -- prefer that for anything new.
    """
    color_color_data = [
        mag_array[0] - mag_array[1], # F1065C - F1140C
        mag_array[1] - mag_array[2], # F1140C - F1550C
        mag_array[0] - mag_array[2]] # F1065C - F1550C
    return color_color_data


def get_full_redshift_mag_loop(file, filter_names=DEFAULT_FILTERS, redshifts=DEFAULT_REDSHIFTS):
    """
    Load `file`, redshift it across `redshifts`, and synth-photometer each step.

    filter_names / redshifts default to the notebook's original hard-coded values
    (the three MIRI coronagraph filters, z = 0.5-2 in 90 steps).

    Returns a list of (z, mag_array, color_color_data) tuples, one per redshift --
    the shape plot_color_color expects.
    """
    wavelengths, fluxes = load_swire_data(file)
    full_mag_data = []
    for z in redshifts:
        red_wavelengths, red_fluxes = redshift_data(wavelengths, fluxes, z)
        mag_data = synth_mags(red_wavelengths, red_fluxes, filter_names)
        # synth_mags keeps the (apparent, error) tuple per filter; drop the error
        un_tupled_data = [item[0] for item in mag_data]
        color_color = get_color_color_data(un_tupled_data)
        full_mag_data.append((z, un_tupled_data, color_color))
    return full_mag_data


def galaxy_color_color_data(file, filter_names=DEFAULT_FILTERS, redshifts=DEFAULT_REDSHIFTS):
    """One call: load a SWIRE template and build every flux/color column across `redshifts`.

    The galaxy-side equivalent of planets'
    update_planet_flux_and_magnitude + add_color_columns pair, collapsed into one
    function since synth_mags is cheap (no ~818 MB HDF5 reads to batch around, unlike
    the planet side). Colors go through the same color_pairs primitive planets use, so
    the "A - B" keys line up and the result is a drop-in argument to
    planets.color_color_matrix.

    Parameters
    ----------
    file : str or Path
        SWIRE template file. Resolve it with
        paths.galaxy_template('swire-library/N6090_template_norm.sed')
        rather than a cwd-relative literal.
    filter_names : sequence of str
        Full species filter names. Defaults to the three MIRI coronagraph filters.
    redshifts : sequence of float
        Redshift grid to evaluate at. Defaults to 90 steps from z=0.5 to 2.

    Returns
    -------
    dict with:
        "file" : str, the input path (for labeling plots/legends)
        "redshift" : ndarray, shape (n_z,)
        "mag_{label}" : ndarray, shape (n_z,), one per filter, e.g. "mag_F1065C"
        "{label_i} - {label_j}" : ndarray, shape (n_z,), one per filter pair
    """
    wavelengths, fluxes = load_swire_data(file)
    labels = [filter_label(name) for name in filter_names]

    mag_cols = {label: [] for label in labels}
    for z in redshifts:
        red_wavelengths, red_fluxes = redshift_data(wavelengths, fluxes, z)
        mags = synth_mags(red_wavelengths, red_fluxes, filter_names)
        # synth_mags keeps the (apparent, error) tuple per filter; drop the error
        for label, (mag, _mag_error) in zip(labels, mags):
            mag_cols[label].append(mag)

    mag_cols = {label: np.asarray(vals, dtype=float) for label, vals in mag_cols.items()}

    out = {"file": str(file), "redshift": np.asarray(redshifts, dtype=float)}
    out.update({f"mag_{label}": arr for label, arr in mag_cols.items()})
    out.update(color_pairs(mag_cols, order=labels))
    return out



def load_k15_data(file): 
    """
    Load the K15 data from a .txt file, which has the following columns
    1. Wavelength (in Microns)
    2. Brightness per unit frequency 
    3. Error bar on the brightness per unit frequency

    The file has three header lines that should be skipped.
    """
    data = np.loadtxt(file, comments='#', skiprows=3, usecols=(0, 1, 2))
    wavelengths = data[:, 0]
    fluxes = data[:, 1]
    errors = data[:, 2]
    return wavelengths, fluxes, errors

def redshift_k15_data(wavelengths, fluxes, z):
    """
    Redshift the K15 data by a given redshift z.

    Parameters:
    - wavelengths: array of wavelengths in microns
    - fluxes: array of fluxes 
    - z: redshift value

    Returns:
    - redshifted_wavelengths: array of redshifted wavelengths in microns
    - redshifted_fluxes: array of redshifted fluxes 
    """
    redshifted_wavelengths = wavelengths * (1 + z)
    redshifted_fluxes = fluxes
    return redshifted_wavelengths, redshifted_fluxes

def translate_k15_L_v_to_f_lambda(wavelengths, fluxes):
    """
    Translate K15 data from L_v to f_lambda.

    Parameters:
    - wavelengths: array of wavelengths in microns
    - fluxes: array of fluxes in L_v

    Returns:
    - translated_fluxes: array of fluxes in f_lambda
    """
    wavelengths_m = wavelengths 
    # Convert L_v to f_lambda using the relation f_lambda = (L_v * c) / (lambda^2)
    c = 2.99792458e14  # Speed of light in microns/s
    translated_fluxes = (fluxes * c) / (wavelengths_m ** 2)
    return translated_fluxes

def _check_filters_fit_template(rest_lo, rest_hi, filter_names, redshifts,
                                template_name):
    """Raise ValueError for any filter a redshifted template cannot cover.

    Shared core of check_filters_fit_k15_templates and
    check_filters_fit_swire_templates -- the redshift geometry and the two
    failure modes are identical for any library, only the loader that produces
    (rest_lo, rest_hi) differs.

    Redshifting multiplies wavelengths by (1 + z), so a template spanning
    [rest_lo, rest_hi] is observed over [rest_lo*(1+z), rest_hi*(1+z)]. Both
    edges move redward with z, so across a redshift grid the tightest blue edge
    is at max(z) and the tightest red edge is at min(z). A filter has to sit
    inside that intersection for every z in the grid, not just some of them.

    Two distinct failures are caught here, and the second is the reason the
    check is strict rather than a warning:

    - No overlap: synth_mags has nothing to integrate and returns NaN. Loud,
      once you find it.
    - Partial overlap: synth_mags integrates a truncated bandpass and returns a
      plausible-looking magnitude that is systematically too faint. Nothing is
      NaN, nothing raises, and every downstream color is quietly wrong. This is
      the failure worth failing fast on.
    """
    z = np.asarray(redshifts, dtype=float)
    z_lo, z_hi = float(z.min()), float(z.max())
    obs_lo = rest_lo * (1.0 + z_hi)
    obs_hi = rest_hi * (1.0 + z_lo)

    for name in filter_names:
        filt_lo, filt_hi = (float(v) for v in ReadFilter(name).wavelength_range())
        if filt_lo >= obs_lo and filt_hi <= obs_hi:
            continue

        if filt_hi <= obs_lo or filt_lo >= obs_hi:
            detail = ("No overlap at all, so synth_mags returns NaN for every "
                      "redshift.")
        else:
            detail = ("Only partial overlap, so synth_mags would integrate a "
                      "truncated bandpass and return a plausible but "
                      "systematically too-faint magnitude -- wrong numbers "
                      "with no error raised.")

        raise ValueError(
            f"{name} ({filt_lo:.3f}-{filt_hi:.3f} um) is not covered by "
            f"{template_name} over z={z_lo:g}-{z_hi:g}: the template spans "
            f"{rest_lo:.3f}-{rest_hi:.1f} um rest-frame, which is observed as "
            f"{obs_lo:.3f}-{obs_hi:.1f} um across that redshift range. {detail} "
            "Drop this filter, or use a template library with wider rest-frame "
            "coverage."
        )


def check_filters_fit_k15_templates(file, filter_names=DEFAULT_FILTERS,
                                    redshifts=DEFAULT_REDSHIFTS):
    """Raise ValueError for any filter the redshifted K15 template cannot cover.

    The galaxy-side counterpart to planets.check_filters_fit_model. Without it a
    coverage mismatch is silent at the point of failure and only surfaces much
    later as something unrecognizable -- adding 2MASS/2MASS.Ks produced
    "LinAlgError: SVD did not converge" from np.linalg.matrix_rank, a hundred
    lines downstream, because the NaN magnitudes rode along in X until something
    finally tried to factor it.

    The K15 libraries all start at 2.0 um rest-frame, which is why nothing
    blueward of roughly 3 um observed can be used with them at z >= 0.5. See
    _check_filters_fit_template for the redshift geometry and the two failure
    modes.
    """
    wavelengths, _fluxes, _errors = load_k15_data(file)
    _check_filters_fit_template(
        float(np.min(wavelengths)), float(np.max(wavelengths)),
        filter_names, redshifts, Path(file).name)


def check_filters_fit_swire_templates(file, filter_names=DEFAULT_FILTERS,
                                      redshifts=DEFAULT_REDSHIFTS):
    """Raise ValueError for any filter the redshifted SWIRE template cannot cover.

    The SWIRE counterpart to check_filters_fit_k15_templates. The Polletta+2007
    templates span ~0.1-6000 um rest-frame, so for the 2MASS/WISE set this never
    actually fires -- but the silent partial-overlap failure it guards against
    is still reachable with a narrow template or an extreme redshift grid, and
    the galaxy_magnitudes_swire cache would otherwise hand back a column of
    plausible, systematically wrong magnitudes with nothing raised.
    """
    wavelengths, _fluxes = load_swire_data(file)
    _check_filters_fit_template(
        float(np.min(wavelengths)), float(np.max(wavelengths)),
        filter_names, redshifts, Path(file).name)


def get_full_redshift_mag_loop_k15(file, filter_names=DEFAULT_FILTERS, redshifts=DEFAULT_REDSHIFTS):
    """
    Load K15 `file`, redshift it across the redshift values, then translate to values species needs, and synthesize per-filter magnitudes.

    Returns a list of (z, mag_array, color_color_data) tuples, one per redshift --
    the shape plot_color_color expects.
    """
    check_filters_fit_k15_templates(file, filter_names, redshifts)

    wavelengths, fluxes, errors = load_k15_data(file)
    # Redshift the data
    full_mag_data = []
    for z in redshifts:
        red_wavelengths, red_fluxes = redshift_k15_data(wavelengths, fluxes, z)
        # Translate to f_lambda
        translated_fluxes = translate_k15_L_v_to_f_lambda(red_wavelengths, red_fluxes)
        mag_data = synth_mags(red_wavelengths, translated_fluxes, filter_names)
        # synth_mags keeps the (apparent, error) tuple per filter; drop the error
        un_tupled_data = [item[0] for item in mag_data]
        color_color = get_color_color_data(un_tupled_data)
        full_mag_data.append((z, un_tupled_data, color_color))
    
    return full_mag_data

def galaxy_color_color_data_k15(file, filter_names=DEFAULT_FILTERS, redshifts=DEFAULT_REDSHIFTS):
    """One call: load a K15 template and build every flux/color column across `redshifts`.

    Parameters
    ----------
    file : str or Path
        K15 template file. Resolve it with
        paths.galaxy_template('K15_templates/MIR_library/MIR0.0.txt')
        rather than a cwd-relative literal.
    filter_names : sequence of str
        Full species filter names. Defaults to the three MIRI coronagraph filters.
    redshifts : sequence of float
        Redshift grid to evaluate at. Defaults to 90 steps from z=0.5 to 2.

    Returns
    -------
    dict with:
        "file" : str, the input path (for labeling plots/legends)
        "redshift" : ndarray, shape (n_z,)
        "mag_{label}" : ndarray, shape (n_z,), one per filter, e.g. "mag_F1065C"
        "{label_i} - {label_j}" : ndarray, shape (n_z,), one per filter pair
    """
    check_filters_fit_k15_templates(file, filter_names, redshifts)

    wavelengths, fluxes, errors = load_k15_data(file)
    labels = [filter_label(name) for name in filter_names]

    mag_cols = {label: [] for label in labels}
    for z in redshifts:
        red_wavelengths, red_fluxes = redshift_k15_data(wavelengths, fluxes, z)
        # Translate to f_lambda
        translated_fluxes = translate_k15_L_v_to_f_lambda(red_wavelengths, red_fluxes)
        mags = synth_mags(red_wavelengths, translated_fluxes, filter_names)
        # synth_mags keeps the (apparent, error) tuple per filter; drop the error
        for label, (mag, _mag_error) in zip(labels, mags):
            mag_cols[label].append(mag)

    result = {"file": str(file), "redshift": np.array(redshifts)}
    for label, mag_list in mag_cols.items():
        result[f"mag_{label}"] = np.array(mag_list)
    for label_i, label_j in combinations(labels, 2):
        result[f"{label_i} - {label_j}"] = result[f"mag_{label_i}"] - result[f"mag_{label_j}"]
    return result