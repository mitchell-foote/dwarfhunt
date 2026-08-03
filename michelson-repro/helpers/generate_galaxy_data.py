"""Galaxy SWIRE-template helpers, moved out of michelson-galaxy-graph.ipynb.

Mirrors generate_planet_list.py's shape: load a template, redshift it, synthesize
per-filter magnitudes, then hand the magnitudes to that module's color_pairs so
planets and galaxies compute colors the same way.
"""

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


def plot_color_color(data, title, ax, color, x_idx, y_idx, x_label, y_label):
    """Scatter one template's color-color track, colored by redshift.

    Accepts either shape this module produces:
    - dict from galaxy_color_color_data -- x/y are looked up directly by
      x_label/y_label (its "A - B" keys), x_idx/y_idx are unused.
    - list of (z, mag_array, color_color_data) tuples from
      get_full_redshift_mag_loop -- x_idx/y_idx index into each row's
      color_color_data (0=F1065C-F1140C, 1=F1140C-F1550C, 2=F1065C-F1550C).

    Feeding the dict shape into the x_idx/y_idx path used to fail with
    "string index out of range": iterating a dict walks its string keys, not
    rows, so item[2] grabbed a single character instead of a color value.
    """
    if isinstance(data, dict):
        z_values = data["redshift"]
        x = data[x_label]
        y = data[y_label]
    else:
        z_values = [item[0] for item in data]
        color_color_data = [item[2] for item in data]
        x = [item[x_idx] for item in color_color_data]
        y = [item[y_idx] for item in color_color_data]

    #ax.scatter(x, y, c=z_values, cmap='viridis', s=30)
    ax.plot(x, y, color=color, alpha=0.3, label=title)
    ax.set_xlabel(x_label); ax.set_ylabel(y_label); ax.legend()


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
