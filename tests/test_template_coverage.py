"""Guards on galaxy-template wavelength coverage.

Two failures hide behind a filter the templates cannot support, and neither
announces itself where it happens:

- No overlap gives NaN magnitudes that ride along in X until something
  unrelated chokes on them. Adding 2MASS/2MASS.Ks surfaced as
  "LinAlgError: SVD did not converge" out of np.linalg.matrix_rank, far from
  the actual cause.
- Partial overlap is worse: synth_mags integrates a truncated bandpass and
  returns a plausible, systematically too-faint magnitude. Nothing is NaN and
  nothing raises, so every downstream color is quietly wrong.

check_filters_fit_k15_templates is the galaxy-side counterpart to
planets.check_filters_fit_model, and these lock in that it catches both.
"""

import numpy as np
import pytest

from dwarfhunt import paths
from dwarfhunt.galaxies import check_filters_fit_k15_templates

TEMPLATE = paths.galaxy_template("K15_templates/MIR_library/MIR0.0.txt")
MIRI = ("JWST/MIRI.F1065C", "JWST/MIRI.F1140C", "JWST/MIRI.F1550C")


def test_mid_ir_filters_are_accepted():
    """The filters the study actually uses must not be blocked."""
    check_filters_fit_k15_templates(
        TEMPLATE, filter_names=MIRI + ("WISE/WISE.W3",),
        redshifts=np.linspace(0.5, 2, 5))


def test_uncovered_filter_is_rejected():
    """Ks sits entirely blueward of the redshifted template -> NaN magnitudes."""
    with pytest.raises(ValueError, match="No overlap"):
        check_filters_fit_k15_templates(
            TEMPLATE, filter_names=("2MASS/2MASS.Ks",),
            redshifts=np.linspace(0.5, 2, 5))


def test_partial_overlap_is_rejected():
    """The dangerous case: a truncated bandpass integrates without complaint.

    At z=0 the template (2.0 um+) covers only the red part of Ks (1.93-2.40 um).
    synth_mags would return a too-faint magnitude rather than NaN, so this must
    raise rather than warn.
    """
    with pytest.raises(ValueError, match="partial overlap"):
        check_filters_fit_k15_templates(
            TEMPLATE, filter_names=("2MASS/2MASS.Ks",), redshifts=np.array([0.0]))


def test_error_names_the_filter_and_the_template():
    """The message has to point at the fix, not just report a failure."""
    with pytest.raises(ValueError) as exc:
        check_filters_fit_k15_templates(
            TEMPLATE, filter_names=("2MASS/2MASS.Ks",),
            redshifts=np.linspace(0.5, 2, 5))
    message = str(exc.value)
    assert "2MASS/2MASS.Ks" in message
    assert "MIR0.0.txt" in message
    assert "um" in message


def test_redshift_range_uses_the_tightest_edge():
    """Coverage must hold for every z in the grid, not merely for some of them.

    F1065C (9.93-11.24 um) is fine at z=0.5 (template starts at 3.0 um observed)
    but not at z=5 (starts at 12.0 um). A grid spanning both must be rejected.
    """
    check_filters_fit_k15_templates(
        TEMPLATE, filter_names=("JWST/MIRI.F1065C",), redshifts=np.array([0.5]))
    with pytest.raises(ValueError):
        check_filters_fit_k15_templates(
            TEMPLATE, filter_names=("JWST/MIRI.F1065C",),
            redshifts=np.array([0.5, 5.0]))
