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
check_filters_fit_swire_templates shares the same core and is exercised the
same way -- its point is that the SWIRE library's ~0.1-6000 um rest-frame span
is what makes the 2MASS/WISE set usable where K15's 2.0 um floor is not.
"""

import numpy as np
import pytest

from dwarfhunt import paths
from dwarfhunt.galaxies import (check_filters_fit_k15_templates,
                                check_filters_fit_swire_templates)

TEMPLATE = paths.galaxy_template("K15_templates/MIR_library/MIR0.0.txt")
MIRI = ("JWST/MIRI.F1065C", "JWST/MIRI.F1140C", "JWST/MIRI.F1550C")

SWIRE = paths.galaxy_template("swire-library/Sc_template_norm.sed")
TWOMASS_WISE = ("2MASS/2MASS.J", "2MASS/2MASS.H", "2MASS/2MASS.Ks",
                "WISE/WISE.W1", "WISE/WISE.W2")


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


# --- SWIRE ------------------------------------------------------------------
# The whole reason for the SWIRE path: 2MASS + WISE through W2 are usable
# across a realistic redshift range, which they are not with K15.


def test_swire_covers_2mass_and_wise_over_a_realistic_redshift_grid():
    check_filters_fit_swire_templates(
        SWIRE, filter_names=TWOMASS_WISE, redshifts=np.linspace(0.0, 3.0, 7))


def test_swire_rejects_no_overlap():
    """Push z high enough that the template's blue edge clears J entirely."""
    with pytest.raises(ValueError, match="No overlap"):
        check_filters_fit_swire_templates(
            SWIRE, filter_names=("2MASS/2MASS.J",),
            redshifts=np.array([10.0, 20.0]))


def test_swire_rejects_partial_overlap():
    """Ks straddles the blue edge at that same z: the silent, too-faint case."""
    with pytest.raises(ValueError, match="partial overlap"):
        check_filters_fit_swire_templates(
            SWIRE, filter_names=("2MASS/2MASS.Ks",),
            redshifts=np.array([10.0, 20.0]))


def test_swire_error_names_the_filter_and_the_template():
    with pytest.raises(ValueError) as exc:
        check_filters_fit_swire_templates(
            SWIRE, filter_names=("2MASS/2MASS.J",),
            redshifts=np.array([10.0, 20.0]))
    message = str(exc.value)
    assert "2MASS/2MASS.J" in message
    assert "Sc_template_norm.sed" in message
    assert "not covered" in message
