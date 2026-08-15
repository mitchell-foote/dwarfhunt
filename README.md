# dwarfhunt

Brown dwarf / galaxy color-color analysis built on [`species`](https://species.readthedocs.io).

The research narrative lives in [Log.md](Log.md). This file is just setup.

## Setup

```bash
conda activate dwarfhunt
pip install -r requirements.txt
pip install -e .
```

The editable install is what makes `import dwarfhunt` work from any directory,
including from inside a notebook in any experiment folder.

## Using it

```python
import dwarfhunt
from dwarfhunt import planets, galaxies, plots, paths

db = dwarfhunt.init()          # attaches to the shared database
```

`init()` generates `species_config.ini` at the repo root on first run, pointing
at `data/species_database.hdf5`. That file holds absolute, machine-specific
paths, so it is gitignored and regenerated per checkout —
`species_config.ini.template` shows its shape.

Pass `force_species_init=True` before any **write** (`add_model`, `add_filter`,
`add_photometry`, `add_companion`), and run those with no other kernels
attached. `SpeciesInit` reopens the database in append mode on every call, so
two kernels writing to one shared database will collide on the HDF5 writer lock.

Galaxy templates resolve through the package rather than a cwd-relative string:

```python
galaxies.galaxy_color_color_data_k15(
    paths.galaxy_template('K15_templates/MIR_library/MIR0.0.txt'))
```

## Layout

```
src/dwarfhunt/     paths, session, planets, galaxies, gmm, plots
data/              shared species database + ~160 GB of model grids  (gitignored)
assets/            galaxy templates (SWIRE, Kirkpatrick+2015)
cache/             derived caches, e.g. the missing-grid-point deny-list
tests/             guards on the shared-config seam
tools/             capture_baseline.py — exact before/after numeric diffs
michelson-repro/   the Michelson reproduction (notebooks)
broadband-filters/ 2MASS / GAIA / WISE profiles, pinned copies from SVO
jwst_filters/      MIRI filter exploration
sprint-week/       early species tutorial work
```

## Working on the data

**`data/` is ~160 GB and the volume runs near capacity. Move it with `mv`,
never `cp`** — a rename is free, a copy will not fit.

The `.tgz` archives in `data/` are not redundant with the extracted `.npy`
directories beside them. `species.add_model` looks for the `.tgz` specifically
and re-downloads it from Leiden when absent, so deleting them costs a ~74 GB
download.

## Tests

```bash
pytest tests/
```

These guard one specific silent failure: if `SPECIES_CONFIG` stops being
exported, species falls back to `$CWD/species_config.ini`, auto-creates an empty
database next to whatever notebook is running, and every number after that is
computed against the wrong grid without raising anything.
