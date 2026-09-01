"""Guards on the missing-spectrum deny-list cache.

The failure this protects against is silent by two hops. `load_deny_list` used
to return the cache file verbatim as soon as it existed, whatever tags had been
asked for -- so once the file held the two Elf Owl tags, asking for any other
tag handed back the Elf Owl entries and no rescan ever happened.

Nothing raises at that point. `generate_planet_arrays` does
`(deny or {}).get(model.model, {})`, finds no entry for its model, and applies
no denial at all; the sample can then land on a grid point species stored as
all zeros, and the only symptom is NaN magnitudes surfacing much later, already
written into the photometry cache.

So the file is checked tag by tag, and the returned dict holds exactly the tags
that were requested -- a missing key is impossible rather than silent.
"""

import json

import pytest

from dwarfhunt import planets


def _entry(name):
    return {"params": [name], "denied_axis_values": {name: [3.0]}, "combos": []}


def test_a_cached_tag_is_reused_without_scanning(tmp_path, monkeypatch):
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("logg")}))

    def explode(*args, **kwargs):
        raise AssertionError("scanned a tag that was already cached")

    monkeypatch.setattr(planets, "_scan_one", explode)
    assert planets.load_deny_list(["tag-a"], path=path) == {"tag-a": _entry("logg")}


def test_a_tag_absent_from_the_file_is_scanned_not_faked(tmp_path, monkeypatch):
    """The actual bug: a warm file must not answer for a tag it never covered."""
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("logg")}))

    monkeypatch.setattr(planets, "_scan_one", lambda tag, db: _entry(f"axis-{tag}"))
    out = planets.load_deny_list(["tag-b"], path=path)

    assert set(out) == {"tag-b"}, "returned a tag that was not requested"
    assert out["tag-b"] == _entry("axis-tag-b")


def test_the_returned_dict_holds_exactly_the_requested_tags(tmp_path, monkeypatch):
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("logg"), "tag-z": _entry("feh")}))

    monkeypatch.setattr(planets, "_scan_one", lambda tag, db: _entry(f"axis-{tag}"))
    out = planets.load_deny_list(["tag-a", "tag-b"], path=path)

    assert set(out) == {"tag-a", "tag-b"}


def test_a_scanned_tag_is_merged_into_the_file_not_replacing_it(tmp_path, monkeypatch):
    """Scanning a tag costs a pass over the whole flux array, so keep the rest."""
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("logg")}))

    monkeypatch.setattr(planets, "_scan_one", lambda tag, db: _entry(f"axis-{tag}"))
    planets.load_deny_list(["tag-b"], path=path)

    assert set(json.loads(path.read_text())) == {"tag-a", "tag-b"}


def test_rebuild_rescans_a_tag_that_was_already_cached(tmp_path, monkeypatch):
    path = tmp_path / "deny.json"
    path.write_text(json.dumps({"tag-a": _entry("stale")}))

    monkeypatch.setattr(planets, "_scan_one", lambda tag, db: _entry("fresh"))
    out = planets.load_deny_list(["tag-a"], path=path, rebuild=True)

    assert out["tag-a"] == _entry("fresh")


def test_unknown_tag_names_the_alternatives(tmp_path):
    """A typo'd tag should say what the database actually holds."""
    with pytest.raises(KeyError, match="Available"):
        planets.scan_missing_grid_points("sonora-not-a-model")
