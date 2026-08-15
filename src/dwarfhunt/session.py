"""One-line bootstrap that points species at the shared database.

Replaces the `SpeciesInit(); db = Database()` preamble. The difference that
matters: after calling `init()` the working directory is irrelevant, so a
notebook can live in any experiment folder and still read the one shared
database.
"""

import os
from pathlib import Path
from typing import Optional, Union

from . import paths


def ensure_config(
    config: Optional[Path] = None,
    *,
    database: Optional[Path] = None,
    data: Optional[Path] = None,
    overwrite: bool = False,
) -> Path:
    """Write the shared species_config.ini if it is missing. Returns its path.

    The file holds absolute, machine-specific paths, so it is gitignored and
    generated per-checkout rather than committed.
    """
    cfg = Path(config) if config is not None else paths.DEFAULT_CONFIG

    if cfg.exists() and not overwrite:
        return cfg

    cfg.parent.mkdir(parents=True, exist_ok=True)
    cfg.write_text(
        paths.render_config(
            database if database is not None else paths.SHARED_DB,
            data if data is not None else paths.SHARED_DATA,
        )
    )
    return cfg


def init(
    config: Optional[Union[str, Path]] = None,
    *,
    force_species_init: bool = False,
    quiet: bool = True,
):
    """Point species at the shared config and return a ready Database.

    Calls `SpeciesInit` only when there is real setup to do -- a missing
    database or `force_species_init`. Otherwise it attaches by exporting
    SPECIES_CONFIG directly, which every species consumer honours identically.

    That distinction is not an optimisation. SpeciesInit reopens the database in
    append mode on every call (it deletes and rewrites the `configuration`
    group), so with one shared database a second notebook kernel calling it
    would collide with the first on the HDF5 writer lock. Attaching takes no
    write lock, and also skips the pypi.org version ping SpeciesInit performs.

    Use `force_species_init=True` before `add_model` / `add_filter` / any other
    write, and run those with no other kernels attached.
    """
    cfg = Path(config) if config is not None else paths.DEFAULT_CONFIG
    ensure_config(cfg)

    # Read the configured database straight from the ini rather than assuming
    # the default, so a caller-supplied config is respected.
    os.environ["SPECIES_CONFIG"] = str(cfg)
    db_path = paths.database_path()

    if force_species_init or not db_path.exists():
        from species import SpeciesInit

        # str(), not Path(). SpeciesInit only exports SPECIES_CONFIG on its
        # isinstance(config_file, str) branch; hand it a Path and it uses the
        # file itself while leaving every downstream class reading the cwd.
        SpeciesInit(config_file=str(cfg))

    if os.environ.get("SPECIES_CONFIG") != str(cfg):
        raise RuntimeError(
            "SPECIES_CONFIG was not exported, so species will read "
            f"{Path.cwd() / 'species_config.ini'} instead of {cfg}. This is the "
            "str-vs-Path trap in SpeciesInit -- every number would be computed "
            "against the wrong database."
        )

    from species.data.database import Database

    database = Database()

    if not quiet:
        print(f"Config:   {cfg}")
        print(f"Database: {database.database}")
        print(f"Data:     {database.data_folder}")

    return database
