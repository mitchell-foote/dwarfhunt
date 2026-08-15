"""Guards on the shared-config mechanism.

The failure these protect against is silent: if SPECIES_CONFIG stops being
exported, species falls back to `os.getcwd()/species_config.ini`, auto-creates
an empty database next to whatever notebook is running, and every subsequent
number is computed against the wrong grid without any error.
"""

import os
from pathlib import Path

import pytest

from dwarfhunt import paths


@pytest.fixture
def clean_env(monkeypatch):
    monkeypatch.delenv("SPECIES_CONFIG", raising=False)


def test_relative_config_values_are_rejected(tmp_path, monkeypatch):
    """A relative path in the ini must fail loudly, not resolve against cwd."""
    cfg = tmp_path / "species_config.ini"
    cfg.write_text(
        "[species]\ndatabase = species_database.hdf5\n"
        "data_folder = ./data/\nvega_mag = 0.03\n"
    )
    monkeypatch.setenv("SPECIES_CONFIG", str(cfg))

    with pytest.raises(RuntimeError, match="relative"):
        paths.database_path()


def test_missing_section_gives_actionable_error(tmp_path, monkeypatch):
    """An absent config yields our message, not species' bare KeyError."""
    monkeypatch.setenv("SPECIES_CONFIG", str(tmp_path / "nope.ini"))

    with pytest.raises(RuntimeError, match="init"):
        paths.database_path()


def test_config_path_follows_species_resolution_order(tmp_path, monkeypatch):
    """config_path() must mirror how species itself picks a config."""
    monkeypatch.delenv("SPECIES_CONFIG", raising=False)
    monkeypatch.chdir(tmp_path)
    assert paths.config_path() == tmp_path / "species_config.ini"

    monkeypatch.setenv("SPECIES_CONFIG", "/somewhere/else.ini")
    assert paths.config_path() == Path("/somewhere/else.ini")


def test_import_is_side_effect_free(tmp_path, monkeypatch):
    """Importing the package must not create files or read the database."""
    monkeypatch.chdir(tmp_path)

    import importlib

    import dwarfhunt

    importlib.reload(dwarfhunt)

    assert list(tmp_path.iterdir()) == []


@pytest.mark.needs_db
@pytest.mark.usefixtures("clean_env")
def test_init_resolves_shared_db_from_foreign_cwd(tmp_path, monkeypatch):
    """The whole point: cwd must not determine which database species reads."""
    if not paths.DEFAULT_CONFIG.exists():
        pytest.skip("shared species_config.ini not generated yet")

    monkeypatch.chdir(tmp_path)

    import dwarfhunt

    database = dwarfhunt.init()

    assert os.environ["SPECIES_CONFIG"] == str(paths.DEFAULT_CONFIG)
    assert Path(database.database) == paths.database_path()
    assert Path(database.database).is_absolute()

    # species must not have auto-created a config beside us
    assert not (tmp_path / "species_config.ini").exists()
