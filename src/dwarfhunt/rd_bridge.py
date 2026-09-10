"""Plumbing between dwarfhunt and Red Dragon (delta).

Red Dragon lives in a separate repo we do not own
(bitbucket.org/wkblack/red-dragon-delta), installed here in editable mode as
`rd_delta`. This module is the seam: format conversion, column selection, and
sentinel scrubbing -- the mechanical parts of getting our arrays into and out
of Red Dragon's file contract.

Deliberately NOT here: the survey depth model, the apparent-magnitude
distribution, and the evidence-ratio scorer. Those are modelling choices, they
belong to the analysis, and they get written in reddragon/ alongside the
notebook that defends them.

Nothing in this module imports species, so it stays cheap.

The colour convention is Red Dragon's and it is not negotiable:

    colors[i] = bands[i] - bands[i+1]      # i.e. -np.diff, bands ordered blue->red

Getting the band order wrong flips the sign of every colour and raises nothing.
"""

from pathlib import Path

import h5py
import numpy as np

# Red Dragon's data contract -- every input file needs exactly these four.
REQUIRED_FIELDS = ("Z", "Z_err", "bands", "bands_err")

# Catalogue "not measured" placeholders seen in the wild. Confirmed present in
# FP_G2MW12_vetted_lgT.h5: every one of its eight bands_err columns maxes at
# exactly 999.
SENTINELS = (999.0, -999.0, 99.0, -99.0)


def bands_to_colors(bands):
    """
    Form colours from magnitudes using Red Dragon's convention.

    Parameters
    ----------
    bands : (N, n_bands) array
        Magnitudes, ordered blue to red.

    Returns
    -------
    (N, n_bands - 1) array of colours, colors[i] = bands[i] - bands[i+1].
    """
    return -np.diff(np.asarray(bands, dtype=float), axis=1)


def scrub_sentinels(arr, sentinels=SENTINELS, tol=1e-6):
    """
    Replace catalogue sentinel values with NaN, so they read as missing.

    Returns a copy; the input is untouched. Red Dragon marginalises over NaN
    colours, but a literal 999 in bands_err is silently accepted as a real
    (enormous) uncertainty, which is worse than missing.
    """
    out = np.array(arr, dtype=float, copy=True)
    for s in sentinels:
        out[np.isclose(out, s, rtol=0, atol=tol)] = np.nan
    return out


def select_bands(path, want, scrub=True):
    """
    Pull a named subset of bands out of a Red Dragon file, in the order asked.

    Reads `label_bands` from the file's attrs and indexes by name rather than
    by position, so a change in the source file's column order cannot silently
    reorder the output.

    Parameters
    ----------
    path : str or Path
        A Red Dragon format HDF5 file carrying a `label_bands` attr.
    want : sequence of str
        Band names to keep, in blue-to-red order, e.g. ("J","H","Ks","W1","W2").
    scrub : bool
        Run scrub_sentinels over bands and bands_err. Default True.

    Returns
    -------
    dict with keys Z, Z_err, bands, bands_err, colors, label_bands, attrs.
    """
    want = list(want)
    with h5py.File(Path(path), "r") as f:
        missing = [k for k in REQUIRED_FIELDS if k not in f]
        if missing:
            raise KeyError(f"{path} is missing required field(s): {missing}")

        if "label_bands" not in f.attrs:
            raise KeyError(
                f"{path} has no label_bands attr, so bands cannot be selected "
                "by name. Index by column only if you have confirmed the order."
            )
        labels = [
            b.decode() if isinstance(b, bytes) else str(b)
            for b in f.attrs["label_bands"]
        ]
        unknown = [w for w in want if w not in labels]
        if unknown:
            raise KeyError(f"band(s) {unknown} not in {labels}")
        idx = [labels.index(w) for w in want]

        Z = f["Z"][:]
        Z_err = f["Z_err"][:]
        bands = f["bands"][:][:, idx]
        bands_err = f["bands_err"][:][:, idx]
        attrs = dict(f.attrs)

    if scrub:
        bands = scrub_sentinels(bands)
        bands_err = scrub_sentinels(bands_err)

    return {
        "Z": Z,
        "Z_err": Z_err,
        "bands": bands,
        "bands_err": bands_err,
        "colors": bands_to_colors(bands),
        "label_bands": want,
        "attrs": attrs,
    }


def write_rd_file(
    path,
    Z,
    Z_err,
    bands,
    bands_err,
    label_primary=None,
    label_bands=None,
    label_colors=None,
    overwrite=False,
):
    """
    Write arrays out in Red Dragon's input format.

    Mechanical transcription of the data contract in red-dragon-delta's
    CLAUDE.md, with the shape checks that catch a transposed array before a
    fit spends an hour on it. Writes to `<path>.tmp` and renames on success,
    matching the never-leave-a-partial-file habit of rd_delta.main.

    Parameters
    ----------
    Z, Z_err : (N,) arrays
        Evolution variable and its uncertainty. Z_err is floored at 1e-16 on
        read by Red Dragon, so an exactly-known Z wants a small positive value
        (1e-4 is the walkthrough's suggestion), not zero.
    bands, bands_err : (N, n_bands) arrays
        Magnitudes ordered blue to red, and their per-band sigmas.
    label_primary, label_bands, label_colors : optional
        Propagated to the fitted model file so plots are self-labelling.

    Raises on zero or non-finite bands_err: passing zero error bars collapses
    extreme deconvolution onto a noiseless manifold, which is the single most
    expensive mistake available here.
    """
    path = Path(path)
    Z = np.asarray(Z, dtype=float)
    Z_err = np.asarray(Z_err, dtype=float)
    bands = np.asarray(bands, dtype=float)
    bands_err = np.asarray(bands_err, dtype=float)

    if Z.ndim != 1:
        raise ValueError(f"Z must be 1-D, got shape {Z.shape}")
    if Z_err.shape != Z.shape:
        raise ValueError(f"Z_err {Z_err.shape} does not match Z {Z.shape}")
    if bands.ndim != 2 or bands.shape[0] != Z.shape[0]:
        raise ValueError(
            f"bands must be (N, n_bands) with N={Z.shape[0]}, got {bands.shape}"
        )
    if bands_err.shape != bands.shape:
        raise ValueError(
            f"bands_err {bands_err.shape} does not match bands {bands.shape}"
        )
    if bands.shape[1] < 2:
        raise ValueError("need at least 2 bands to form a colour")

    # The zero-error trap. NaN is fine -- that means missing, and Red Dragon
    # marginalises over it. Zero is not: it claims a perfect measurement.
    finite = np.isfinite(bands_err)
    if np.any(bands_err[finite] <= 0):
        n = int(np.sum(bands_err[finite] <= 0))
        raise ValueError(
            f"{n} bands_err entries are <= 0. Red Dragon's likelihood uses "
            "Sigma + Delta_j; with Delta_j = 0 the fit drives components "
            "infinitely narrow against a noiseless manifold. Use a real depth "
            "model, or NaN for genuinely missing bands."
        )

    if label_bands is not None and len(label_bands) != bands.shape[1]:
        raise ValueError(
            f"label_bands has {len(label_bands)} entries but bands has "
            f"{bands.shape[1]} columns"
        )
    if label_colors is not None and len(label_colors) != bands.shape[1] - 1:
        raise ValueError(
            f"label_colors has {len(label_colors)} entries but "
            f"{bands.shape[1]} bands give {bands.shape[1] - 1} colours"
        )

    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists; pass overwrite=True to replace")

    tmp = path.with_suffix(path.suffix + ".tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(tmp, "w") as ds:
        ds.create_dataset("Z", data=Z)
        ds.create_dataset("Z_err", data=Z_err)
        ds.create_dataset("bands", data=bands)
        ds.create_dataset("bands_err", data=bands_err)
        if label_primary is not None:
            ds.attrs["label_primary"] = label_primary
        if label_bands is not None:
            ds.attrs["label_bands"] = list(label_bands)
        if label_colors is not None:
            ds.attrs["label_colors"] = list(label_colors)
    tmp.rename(path)
    return path


def colors_and_covars(bands, bands_err):
    """
    Colours plus their full tridiagonal covariance, ready for get_posterior_Z.

    Adjacent colours share a band, so the colour covariance carries
    -sigma_band^2 off-diagonals. Dropping them is wrong, and rd_delta's fitting
    code depends on their being right -- this just routes to the repo's own
    builder so there is one implementation.

    Returns
    -------
    colors : (N, D) array
    covars : (N, D, D) array
    """
    from rd_delta.utils import band_err_to_covar

    return bands_to_colors(bands), band_err_to_covar(np.asarray(bands_err, float))
