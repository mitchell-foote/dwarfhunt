# dwarfhunt

Brown dwarf / galaxy color-color analysis built on [`species`](https://species.readthedocs.io).

## Start here

**[`pop-separation/pop-separation.ipynb`](pop-separation/pop-separation.ipynb)**
is the entry point for the project. It separates Sonora Bobcat brown dwarf models from Kirkpatrick+2015 galaxy templates in JWST/MIRI color space using a Gaussian-mixture classifier, and walks through, in order:

1. Building the dwarf and galaxy color samples
2. Choosing the number of GMM components by BIC
3. The confusion matrix and decision regions for each color pair
4. Which dwarfs (by T_eff) and galaxies (by redshift / AGN fraction) are misclassified
5. Whether each extra color improves balanced accuracy, paired across 20 splits
6. A filter-subset search, scored once on a sealed holdout
7. Injected photometric noise, and the scatter at which separation breaks down

The notebook is saved with its outputs, so you can read it without running anything. Re-running it needs the Sonora Bobcat grid in the species database (see Setup below).

The research log lives in [Log.md](Log.md).

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
from dwarfhunt import dwarfs, galaxies, plots, paths

db = dwarfhunt.init()          # attaches to the shared database
```

`init()` generates `species_config.ini` at the repo root on first run, pointing at `data/species_database.hdf5`. That file holds absolute, machine-specific paths, so it is gitignored and regenerated per checkout — `species_config.ini.template` shows its shape.

Pass `force_species_init=True` before any **write** (`add_model`, `add_filter`, `add_photometry`, `add_companion`), and run those with no other kernels attached. 

`SpeciesInit` reopens the database in append mode on every call, so two kernels writing to one shared database will collide on the HDF5 writer lock.

Galaxy templates resolve through the package rather than a cwd-relative string:

```python
galaxies.galaxy_color_color_data_k15(
    paths.galaxy_template('K15_templates/MIR_library/MIR0.0.txt'))
```

## Layout

```
src/dwarfhunt/     paths, session, dwarfs, galaxies, gmm, plots
data/              shared species database + ~160 GB of model grids  (gitignored)
assets/            galaxy templates (SWIRE, Kirkpatrick+2015)
cache/             derived caches, e.g. the missing-grid-point deny-list
tests/             guards on the shared-config seam
tools/             capture_baseline.py — exact before/after numeric diffs
pop-separation/    ** entry point ** dwarf/galaxy separation in MIRI colors
2mass-wise-separation/  follow-on: 2MASS + WISE filters, Elf Owl / Diamondback + SWIRE
michelson-repro/   the Michelson reproduction (notebooks)
broadband-filters/ 2MASS / GAIA / WISE profiles, pinned copies from SVO
jwst_filters/      MIRI filter exploration
sprint-week/       early species tutorial work
```

## Working on the data

**`data/` directory can get ~160 GB with all downloads, so if you want to run the models locally, make sure you have space on your disk. 

The `.tgz` archives in `data/` are not redundant with the extracted `.npy` directories beside them. `species.add_model` looks for the `.tgz` specifically and re-downloads it from Leiden when absent, so deleting them costs a ~74 GB download.

## Tests

```bash
pytest tests/
```
