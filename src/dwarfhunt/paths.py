"""Where everything lives.

This module is the bridge between two path systems that must never disagree:
the one `species` uses internally, and the one our own direct-HDF5 readers use
when they bypass species and index the grid arrays themselves.

The rule here is that the ini file is the authority, not this module. species
re-reads `species_config.ini` independently in 20+ classes and there is no way
to make it consult us instead, so we read the same file it does rather than
keeping a second copy of the paths. `ensure_config` writes that file on first
run; after that it is just read.

Only stdlib is imported, deliberately: the direct-HDF5 readers need a database
path without dragging in species or triggering SpeciesInit.
"""

import os
from configparser import ConfigParser
from pathlib import Path

# src/dwarfhunt/paths.py -> src/dwarfhunt -> src -> repo root. Valid because the
# package is installed editable, so the source stays in the repo.
REPO_ROOT = Path(os.environ.get("DWARFHUNT_ROOT", Path(__file__).resolve().parents[2]))

DEFAULT_CONFIG = REPO_ROOT / "species_config.ini"
SHARED_DATA = REPO_ROOT / "data"
# Inside the data folder rather than beside it. Safe because nothing in species
# enumerates data_folder wholesale -- it only ever joins a known model tag onto
# it -- so the database sitting there is never mistaken for grid content.
SHARED_DB = SHARED_DATA / "species_database.hdf5"
GALAXY_TEMPLATES = REPO_ROOT / "assets" / "galaxy-templates"
CACHE_DIR = REPO_ROOT / "cache"

VEGA_MAG = "0.03"


def config_path() -> Path:
    """The ini species itself will read: SPECIES_CONFIG if set, else cwd.

    Mirrors the resolution order in species.data.database.Database.__init__ and
    its ~20 siblings, so that asking this module which config is live gives the
    same answer species would give.
    """
    return Path(os.environ.get("SPECIES_CONFIG", Path.cwd() / "species_config.ini"))


def _config_value(key: str) -> Path:
    cfg = config_path()
    parser = ConfigParser()
    parser.read(cfg)

    try:
        raw = parser["species"][key]
    except KeyError:
        raise RuntimeError(
            f"{cfg} has no [species] {key}. Call dwarfhunt.init() before anything "
            "that reads the database."
        ) from None

    value = Path(raw)

    if not value.is_absolute():
        # species applies no normalisation anywhere -- the raw string goes
        # straight to h5py.File -- so a relative value resolves against whatever
        # cwd the caller happens to have. That is the mechanism behind the
        # silently-reading-the-wrong-database class of bug, so refuse it here
        # rather than let this module and species resolve it differently.
        raise RuntimeError(
            f"{cfg} sets {key}={raw!r}, which is relative. Use an absolute path: "
            "species resolves relative values against the current working "
            "directory, not against the config file's location."
        )

    return value


def database_path() -> Path:
    """Absolute path to the species database that is currently configured."""
    return _config_value("database")


def data_folder() -> Path:
    """Absolute path to the species download cache that is currently configured."""
    return _config_value("data_folder")


def galaxy_template(relative: str) -> Path:
    """Resolve a galaxy template path, e.g. 'K15_templates/MIR_library/MIR0.0.txt'."""
    return GALAXY_TEMPLATES / relative


def cache_file(name: str) -> Path:
    """Resolve a path under the repo cache directory, creating the directory."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / name


def render_config(database: Path, data: Path) -> str:
    """The text of a species_config.ini pointing at `database` and `data`."""
    return (
        "[species]\n"
        f"database = {database}\n"
        f"data_folder = {data}\n"
        f"vega_mag = {VEGA_MAG}\n"
    )
