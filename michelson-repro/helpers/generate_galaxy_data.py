"""Galaxy template helpers (SWIRE and Kirkpatrick+2015), moved out of
michelson-galaxy-graph.ipynb.

Mirrors generate_planet_list.py's shape: load a template, redshift it, synthesize
per-filter magnitudes, then hand the magnitudes to that module's color_pairs so
planets and galaxies compute colors the same way.

Data only -- the drawing lives in color_color_plots.py.
"""

from itertools import combinations

import numpy as np
from species.phot.syn_phot import SyntheticPhotometry

from .generate_planet_list import color_pairs, filter_label

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

    The galaxy-side equivalent of generate_planet_list's
    update_planet_flux_and_magnitude + add_color_columns pair, collapsed into one
    function since synth_mags is cheap (no ~818 MB HDF5 reads to batch around, unlike
    the planet side). Colors go through the same color_pairs primitive planets use, so
    the "A - B" keys line up and the result is a drop-in argument to
    generate_planet_list.color_color_matrix.

    Parameters
    ----------
    file : str or Path
        SWIRE template file, e.g. "./galaxy-data/swire-library/N6090_template_norm.sed".
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

def synth_mags_k15(wavelength, flux, filter_names):
    """
    Calculate synthetic magnitudes for a given list of filters using K15 data.

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

def get_full_redshift_mag_loop_k15(file, filter_names=DEFAULT_FILTERS, redshifts=DEFAULT_REDSHIFTS):
    """
    Load K15 `file`, redshift it across the redshift values, then translate to values species needs, and synthesize per-filter magnitudes.

    Returns a list of (z, mag_array, color_color_data) tuples, one per redshift --
    the shape plot_color_color expects.
    """
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
        K15 template file, e.g. "./galaxy-data/K15_templates/MIR_library/MIRO0.0.txt".
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