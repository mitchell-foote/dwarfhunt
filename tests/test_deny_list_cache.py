"""Guards on the missing-spectrum deny-list cache.

Two failures hide behind this cache, and both are silent unless caught here.

`load_deny_list` used to return the cache file verbatim as soon as it existed,
whatever tags had been asked for -- so once the file held the two Elf Owl
tags, asking for any other tag handed back the Elf Owl entries and no rescan
ever happened. Nothing raises at that point: `generate_planet_arrays` does
`(deny or {}).get(model.model, {})`, finds no entry for its model, and applies
no denial at all; the sample can then land on a grid point species stored as
all zeros, and the only symptom is NaN magnitudes surfacing much later,
already written into the photometry cache. So the file is checked tag by tag,
and the returned dict holds exactly the tags that were requested -- a missing
key is impossible rather than silent.

The second failure is the same shape one layer up: a cached entry can go
stale. Which grid points are missing does not change when a model is re-added
with a different wavel_range/wavel_sampling (species resamples an existing
point, it doesn't add or remove one), but it can change for other reasons -- a
teff_range subset, a re-extraction, an updated species version reusing the
same tag -- and a warm cache has no way to notice on its own. So every cached
tag's `grid_shape` (the stored flux array's shape, cheap to read) is compared
against the live database before it is trusted; a mismatch is treated exactly
like the tag never having been cached.
"""

import json

import pytest

from dwarfhunt import planets


def _entry(name, grid_shape=(1, 1)):
    return {"params": [name], "denied_axis_values": {name: [3.0]}, "combos": [],
            "grid_shape": list(grid_shape)}


def _fixed_shape(monkeypatch, shape=(1, 1)):
    """Stand in for _grid_shape so these tests never touch a real database."""
    monkeypatch.setattr(planets, "_grid_shape", lambda tag, db: list(shape))


def test_a_cached_tag_is_reused_without_scanning(tmp_path, monkeypatch):
    _fixed_shape(monkeypatch)
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("logg")}))

    def explode(*args, **kwargs):
        raise AssertionError("scanned a tag that was already cached")

    monkeypatch.setattr(planets, "_scan_one", explode)
    assert planets.load_deny_list(["tag-a"], path=path) == {"tag-a": _entry("logg")}


def test_a_tag_absent_from_the_file_is_scanned_not_faked(tmp_path, monkeypatch):
    """The actual bug: a warm file must not answer for a tag it never covered."""
    _fixed_shape(monkeypatch)
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("logg")}))

    monkeypatch.setattr(planets, "_scan_one", lambda tag, db: _entry(f"axis-{tag}"))
    out = planets.load_deny_list(["tag-b"], path=path)

    assert set(out) == {"tag-b"}, "returned a tag that was not requested"
    assert out["tag-b"] == _entry("axis-tag-b")


def test_the_returned_dict_holds_exactly_the_requested_tags(tmp_path, monkeypatch):
    _fixed_shape(monkeypatch)
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("logg"), "tag-z": _entry("feh")}))

    monkeypatch.setattr(planets, "_scan_one", lambda tag, db: _entry(f"axis-{tag}"))
    out = planets.load_deny_list(["tag-a", "tag-b"], path=path)

    assert set(out) == {"tag-a", "tag-b"}


def test_a_scanned_tag_is_merged_into_the_file_not_replacing_it(tmp_path, monkeypatch):
    """Scanning a tag costs a pass over the whole flux array, so keep the rest."""
    _fixed_shape(monkeypatch)
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("logg")}))

    monkeypatch.setattr(planets, "_scan_one", lambda tag, db: _entry(f"axis-{tag}"))
    planets.load_deny_list(["tag-b"], path=path)

    assert set(json.loads(path.read_text())) == {"tag-a", "tag-b"}


def test_rebuild_rescans_a_tag_that_was_already_cached(tmp_path, monkeypatch):
    _fixed_shape(monkeypatch)
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("stale")}))

    monkeypatch.setattr(planets, "_scan_one", lambda tag, db: _entry("fresh"))
    out = planets.load_deny_list(["tag-a"], path=path, rebuild=True)

    assert out["tag-a"] == _entry("fresh")


def test_unknown_tag_names_the_alternatives(tmp_path):
    """A typo'd tag should say what the database actually holds."""
    with pytest.raises(KeyError, match="Available"):
        planets.scan_missing_grid_points("sonora-not-a-model")


# --- staleness: a cached entry that no longer matches the live grid --------


def test_a_shape_mismatch_forces_a_rescan(tmp_path, monkeypatch):
    """The new failure this guards: re-adding a model must not leave a stale
    deny-list in play just because its tag was already in the cache file."""
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("logg", grid_shape=(1, 1, 100))}))

    monkeypatch.setattr(planets, "_grid_shape", lambda tag, db: [1, 1, 5222])
    monkeypatch.setattr(planets, "_scan_one",
                        lambda tag, db: _entry("logg", grid_shape=(1, 1, 5222)))

    out = planets.load_deny_list(["tag-a"], path=path)
    assert out["tag-a"] == _entry("logg", grid_shape=(1, 1, 5222))


def test_a_shape_match_is_reused_without_scanning(tmp_path, monkeypatch):
    """The companion case: a live grid that still matches must NOT be rescanned
    on every call -- scanning costs a pass over the whole flux array."""
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("logg", grid_shape=(1, 1, 5222))}))

    monkeypatch.setattr(planets, "_grid_shape", lambda tag, db: [1, 1, 5222])
    monkeypatch.setattr(
        planets, "_scan_one",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("scanned a tag whose shape still matched")))

    out = planets.load_deny_list(["tag-a"], path=path)
    assert out["tag-a"] == _entry("logg", grid_shape=(1, 1, 5222))


def test_an_entry_written_before_this_check_existed_is_treated_as_stale(tmp_path, monkeypatch):
    """A cache file from before grid_shape was tracked has no field to compare --
    that must fail safe (rescan), not read as `None == None` and pass."""
    path = tmp_path / "deny.json"
    old_entry = {"params": ["logg"], "denied_axis_values": {"logg": [3.0]}, "combos": []}
    path.write_text(json.dumps({"tag-a": old_entry}))

    monkeypatch.setattr(planets, "_grid_shape", lambda tag, db: [1, 1, 5222])
    monkeypatch.setattr(planets, "_scan_one",
                        lambda tag, db: _entry("logg", grid_shape=(1, 1, 5222)))

    out = planets.load_deny_list(["tag-a"], path=path)
    assert out["tag-a"]["grid_shape"] == [1, 1, 5222]
