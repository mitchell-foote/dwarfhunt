"""Guards on the dwarfhunt -> Red Dragon format seam.

Two failure modes here are silent, and both are why this file exists:

  * a wrong band order flips the sign of every colour and raises nothing, so
    select_bands indexes by name and this asserts it does;
  * zero bands_err collapses extreme deconvolution onto a noiseless manifold,
    producing a plausible-looking model that is wrong, so write_rd_file
    refuses and this asserts it refuses.

No species database needed; these run anywhere.
"""

import numpy as np
import pytest

from dwarfhunt import rd_bridge as rb

h5py = pytest.importorskip("h5py")


@pytest.fixture
def eight_band(tmp_path):
    """A miniature stand-in for FP_G2MW12_vetted_lgT.h5."""
    rng = np.random.default_rng(0)
    n = 12
    path = tmp_path / "eight.h5"
    bands = rng.uniform(12.0, 18.0, (n, 8))
    bands_err = np.full((n, 8), 0.03)
    bands_err[0, 5] = 999.0  # sentinel with a valid magnitude beside it
    bands[1, 3] = np.nan  # genuinely missing magnitude
    with h5py.File(path, "w") as f:
        f.create_dataset("Z", data=rng.uniform(2.6, 3.5, n))
        f.create_dataset("Z_err", data=np.full(n, 0.01))
        f.create_dataset("bands", data=bands)
        f.create_dataset("bands_err", data=bands_err)
        f.attrs["label_bands"] = ["BP", "G", "RP", "J", "H", "Ks", "W1", "W2"]
    return path, bands


def test_colors_follow_the_minus_diff_convention():
    bands = np.array([[1.0, 2.0, 4.0]])
    assert np.allclose(rb.bands_to_colors(bands), [[-1.0, -2.0]])


def test_select_bands_indexes_by_name_not_position(eight_band):
    path, bands = eight_band
    d = rb.select_bands(path, ("J", "H", "Ks", "W1", "W2"))
    assert d["bands"].shape == (12, 5)
    assert d["colors"].shape == (12, 4)
    # columns 3..7 of the source, selected by label
    np.testing.assert_allclose(d["bands"][:, 0], bands[:, 3], equal_nan=True)
    np.testing.assert_allclose(d["bands"][:, -1], bands[:, 7], equal_nan=True)


def test_select_bands_honours_requested_order(eight_band):
    path, bands = eight_band
    d = rb.select_bands(path, ("W2", "J"))
    np.testing.assert_allclose(d["bands"][:, 0], bands[:, 7], equal_nan=True)
    np.testing.assert_allclose(d["bands"][:, 1], bands[:, 3], equal_nan=True)


def test_select_bands_rejects_unknown_band(eight_band):
    path, _ = eight_band
    with pytest.raises(KeyError, match="Y"):
        rb.select_bands(path, ("J", "Y"))


def test_sentinels_become_nan(eight_band):
    path, _ = eight_band
    d = rb.select_bands(path, ("J", "H", "Ks", "W1", "W2"))
    assert np.isnan(d["bands_err"][0, 2])  # the 999 in Ks
    assert np.isfinite(d["bands_err"][1, 2])


def test_scrub_does_not_mutate_input():
    arr = np.array([[0.03, 999.0]])
    out = rb.scrub_sentinels(arr)
    assert np.isnan(out[0, 1])
    assert arr[0, 1] == 999.0


def test_write_rejects_zero_errors(tmp_path):
    n = 5
    with pytest.raises(ValueError, match="bands_err"):
        rb.write_rd_file(
            tmp_path / "z.h5",
            np.linspace(0, 1, n),
            np.full(n, 1e-4),
            np.zeros((n, 3)),
            np.zeros((n, 3)),
        )


def test_write_allows_nan_errors_as_missing(tmp_path):
    n = 5
    err = np.full((n, 3), 0.05)
    err[0, 1] = np.nan
    p = rb.write_rd_file(
        tmp_path / "nan.h5",
        np.linspace(0, 1, n),
        np.full(n, 1e-4),
        np.full((n, 3), 15.0),
        err,
    )
    assert p.exists()


def test_write_rejects_mismatched_shapes(tmp_path):
    with pytest.raises(ValueError, match="bands_err"):
        rb.write_rd_file(
            tmp_path / "m.h5",
            np.zeros(5),
            np.full(5, 1e-4),
            np.full((5, 3), 15.0),
            np.full((5, 4), 0.05),
        )


def test_write_rejects_label_count_mismatch(tmp_path):
    with pytest.raises(ValueError, match="label_colors"):
        rb.write_rd_file(
            tmp_path / "l.h5",
            np.zeros(5),
            np.full(5, 1e-4),
            np.full((5, 3), 15.0),
            np.full((5, 3), 0.05),
            label_colors=["a", "b", "c"],  # 3 bands give 2 colours
        )


def test_write_round_trips_the_contract(tmp_path):
    n = 7
    rng = np.random.default_rng(1)
    Z = rng.uniform(0, 8, n)
    bands = rng.uniform(14, 18, (n, 5))
    err = np.full((n, 5), 0.05)
    p = rb.write_rd_file(
        tmp_path / "gal.h5",
        Z,
        np.full(n, 1e-4),
        bands,
        err,
        label_primary=r"$z$",
        label_bands=["J", "H", "Ks", "W1", "W2"],
        label_colors=["a", "b", "c", "d"],
    )
    with h5py.File(p, "r") as f:
        assert set(rb.REQUIRED_FIELDS) <= set(f.keys())
        np.testing.assert_allclose(f["Z"][:], Z)
        np.testing.assert_allclose(f["bands"][:], bands)
        assert list(f.attrs["label_bands"]) == ["J", "H", "Ks", "W1", "W2"]
    assert not p.with_suffix(p.suffix + ".tmp").exists()


def test_write_refuses_to_clobber(tmp_path):
    n = 3
    args = (np.zeros(n), np.full(n, 1e-4), np.full((n, 3), 15.0), np.full((n, 3), 0.05))
    p = tmp_path / "once.h5"
    rb.write_rd_file(p, *args)
    with pytest.raises(FileExistsError):
        rb.write_rd_file(p, *args)
    rb.write_rd_file(p, *args, overwrite=True)


def test_covars_are_tridiagonal_symmetric_and_psd():
    rng = np.random.default_rng(2)
    err = rng.uniform(0.01, 0.1, (6, 5))
    colors, covars = rb.colors_and_covars(rng.uniform(14, 18, (6, 5)), err)
    assert colors.shape == (6, 4)
    assert covars.shape == (6, 4, 4)
    assert np.allclose(covars, np.swapaxes(covars, 1, 2))
    assert (np.linalg.eigvalsh(covars) > 0).all()
    # adjacent colours share a band: off-diagonals are -sigma_band^2
    np.testing.assert_allclose(covars[:, 0, 1], -err[:, 1] ** 2)
    # non-adjacent colours are uncorrelated
    assert np.allclose(covars[:, 0, 2], 0.0)
