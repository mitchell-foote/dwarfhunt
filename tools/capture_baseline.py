"""Freeze the numbers the restructure must not change, and diff them later.

Usage
-----
    # before touching anything
    python tools/capture_baseline.py capture -o baseline_before.npz

    # after each migration phase
    python tools/capture_baseline.py capture -o baseline_after.npz
    python tools/capture_baseline.py compare baseline_before.npz baseline_after.npz

`compare` exits non-zero if anything moved, so it can gate a phase.

Why a script and not a notebook: this has to produce identical output on both
sides of the move, and a notebook carries kernel state and a working directory
that will not be the same before and after.

Why exact comparison and not a tolerance: nothing in the restructure is
*supposed* to change a bit. Moving files and renaming modules cannot legitimately
perturb a float. A tolerance would hide exactly the drift this is hunting for.

The one phase that DOES change numbers on purpose is the bug-fix phase (n-dim
grid indexing and the searchsorted fix). Re-capture a fresh `before` immediately
prior to that phase so the intentional deltas are isolated from the move.
"""

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]

# Fixed everywhere. Changing these invalidates an existing baseline.
SEED = 0
N_PLANETS = 50
RADIUS_RANGE = (0.6, 1.3)
DISTANCE = 10

# One 3-axis grid and one 5-axis grid. Both are deliberate: the direct-HDF5
# reader indexes grid axes positionally, so a 3-param grid and a 5-param grid
# exercise different code paths through it.
GRIDS = ("sonora-bobcat", "sonora-elfowl-t")

MIRI = ("JWST/MIRI.F1065C", "JWST/MIRI.F1140C", "JWST/MIRI.F1550C")


def load_helpers():
    """Import the analysis code from wherever it currently lives.

    The same script has to run before the restructure (code in
    michelson-repro/helpers/) and after it (code in src/dwarfhunt/). Trying the
    package first and falling back keeps one script valid on both sides, which
    is what makes the before/after diff meaningful.
    """
    try:
        from dwarfhunt import galaxies, gmm, planets, plots  # noqa: F401

        return {"planets": planets, "grids": planets, "colors": planets}
    except ImportError:
        sys.path.insert(0, str(REPO_ROOT / "michelson-repro"))
        import helpers.generate_planet_list as legacy

        # Pre-split, one module plays all three roles.
        return {"planets": legacy, "grids": legacy, "colors": legacy}


def make_rng():
    """A fresh seeded Generator, so every capture_* starts from the same stream."""
    return np.random.default_rng(SEED)


# --------------------------------------------------------------------------
# Worked examples -- these two are complete
# --------------------------------------------------------------------------


def capture_sampling(h, db):
    """The planet draw itself, per grid.

    Cheap, and the foundation under everything else: if the sampler drifts,
    every downstream number drifts with it and this pins down why.
    """
    from species.read.read_model import ReadModel

    out = {}

    for tag in GRIDS:
        model = ReadModel(tag)
        arrays = h["planets"].generate_planet_arrays(
            model,
            radius_range=RADIUS_RANGE,
            distance=DISTANCE,
            num_samples=N_PLANETS,
            rng=make_rng(),
        )
        for key, values in arrays.items():
            out[f"sampling/{tag}/{key}"] = np.asarray(values, dtype=float)

    return out


def capture_spectra_agreement(h, db):
    """Direct-HDF5 spectra vs the species path, for both grid shapes.

    This is the existing cross-check in the codebase promoted to a gate. It is
    the regression test for any change to the direct reader's axis indexing.
    """
    from species.read.read_model import ReadModel

    out = {}

    for tag in GRIDS:
        model = ReadModel(tag)
        arrays = h["planets"].generate_planet_arrays(
            model,
            radius_range=RADIUS_RANGE,
            distance=DISTANCE,
            num_samples=5,
            rng=make_rng(),
        )

        lo, hi = h["grids"].model_wavel_range(tag)
        window = (lo, min(hi, 14.9))

        try:
            wl, flux = h["grids"].get_planet_spectra(tag, arrays, wavel_range=window)
            out[f"spectra/{tag}/wavelength"] = wl
            out[f"spectra/{tag}/flux"] = flux
        except (IndexError, ValueError) as exc:
            # Expected on 3-axis grids until the direct reader is generalised.
            # Recorded rather than raised so the baseline still captures.
            print(f"  note: get_planet_spectra failed on {tag}: {exc}")

    return out


# --------------------------------------------------------------------------
# TODO which further numbers are load-bearing
# --------------------------------------------------------------------------


def capture_photometry(h, db):
    """TODO: fluxes, absolute magnitudes and colours per grid.

    Sketch:
        arrays = h["planets"].generate_planet_arrays(model, ..., rng=make_rng())
        arrays = h["planets"].update_planet_flux_and_magnitude(model, arrays, MIRI)
        table  = h["colors"].add_color_columns(arrays, MIRI)
        names, matrix = h["colors"].color_color_matrix(table)
        out["colors/<tag>/matrix"] = matrix

    Note F1550C reaches 16.2 um, past the elfowl grid edge -- decide whether
    this covers bobcat only, or both grids with a reduced filter set.
    """
    return {}


def capture_extrapolation(h, db):
    """TODO: the power-law residual tables.

    Sketch:
        result = h["grids"].extrapolation_residuals(
            tag, arrays, fit_range=(12.0, 14.0), test_range=(14.0, 14.9))
        out[f"residual/{tag}"] = result["residual"]

    These are the numbers in Log.md's 2026-07-30 and 07-31 tables, so drift here
    would contradict written conclusions.
    """
    return {}


def capture_galaxies(h, db):
    """TODO: galaxy template colours.

    Sketch, using the shared template root so it is move-invariant:
        from dwarfhunt import paths
        data = galaxies.galaxy_color_color_data_k15(
            paths.galaxy_template('K15_templates/MIR_library/MIR0.0.txt'))

    Worth including if you want the galaxies dedup gated -- that change touches
    the code path these run through.
    """
    return {}


def capture_gmm(h, db):
    """TODO: the classifier numbers.

    Sketch:
        from dwarfhunt.gmm import GMMClassifier, balanced_accuracy
        X, y = <the color matrix and labels>
        for k in range(2, 9):
            clf = GMMClassifier(n_components=k, n_init=10, random_state=SEED).fit(X, y)
            bic.append(clf.bic(X)); acc.append(clf.score(X, y))

    Needs a fixed train/test split -- pass random_state=SEED to
    train_test_split or the accuracy will not reproduce.
    """
    return {}


CAPTURES = (
    capture_sampling,
    capture_spectra_agreement,
    capture_photometry,
    capture_extrapolation,
    capture_galaxies,
    capture_gmm,
)


# --------------------------------------------------------------------------
# Mechanics -- no need to touch below here
# --------------------------------------------------------------------------


def capture(output: Path) -> None:
    import dwarfhunt

    db = dwarfhunt.init()
    print(f"Database: {db.database}\n")

    h = load_helpers()
    arrays: dict[str, np.ndarray] = {}

    for fn in CAPTURES:
        name = fn.__name__.replace("capture_", "")
        produced = fn(h, db)
        if not produced:
            print(f"[skip] {name} (not filled in)")
            continue
        overlap = set(produced) & set(arrays)
        if overlap:
            raise KeyError(f"{name} reuses existing keys: {sorted(overlap)}")
        arrays.update(produced)
        print(f"[ok]   {name}: {len(produced)} arrays")

    if not arrays:
        raise SystemExit(
            "\nNothing captured. Fill in at least one capture_* function before "
            "using this as a gate -- an empty baseline passes every comparison."
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **arrays)
    print(f"\nWrote {len(arrays)} arrays to {output}")

    for key in sorted(arrays):
        print(f"  {digest(arrays[key])}  {key}")


def digest(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()[:12]


def identical(a: np.ndarray, b: np.ndarray) -> bool:
    if a.shape != b.shape or a.dtype != b.dtype:
        return False
    if np.issubdtype(a.dtype, np.floating):
        # NaN marks a missing spectrum in the flux/mag columns, and a plain
        # equality check would call every NaN a difference.
        return np.array_equal(a, b, equal_nan=True)
    return np.array_equal(a, b)


def compare(before: Path, after: Path) -> int:
    left, right = np.load(before), np.load(after)
    lkeys, rkeys = set(left.files), set(right.files)

    problems = []

    for key in sorted(lkeys - rkeys):
        problems.append(f"MISSING  {key}  (in {before.name}, absent from {after.name})")
    for key in sorted(rkeys - lkeys):
        problems.append(f"ADDED    {key}  (absent from {before.name})")

    for key in sorted(lkeys & rkeys):
        a, b = left[key], right[key]
        if identical(a, b):
            continue
        detail = f"shape {a.shape}->{b.shape}" if a.shape != b.shape else ""
        if not detail and np.issubdtype(a.dtype, np.floating):
            with np.errstate(invalid="ignore"):
                worst = np.nanmax(np.abs(a - b)) if a.size else float("nan")
            detail = f"max abs diff {worst:.3e}"
        problems.append(f"CHANGED  {key}  {detail}")

    if problems:
        print(f"{len(problems)} difference(s):\n")
        for line in problems:
            print(f"  {line}")
        print("\nThese should be bit-identical. Investigate before continuing.")
        return 1

    print(f"All {len(lkeys)} arrays identical.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    cap = sub.add_parser("capture", help="write a baseline archive")
    cap.add_argument("-o", "--output", type=Path, default=Path("baseline.npz"))

    cmp_ = sub.add_parser("compare", help="diff two baseline archives")
    cmp_.add_argument("before", type=Path)
    cmp_.add_argument("after", type=Path)

    args = parser.parse_args()

    if args.command == "capture":
        capture(args.output)
        return 0
    return compare(args.before, args.after)


if __name__ == "__main__":
    raise SystemExit(main())
