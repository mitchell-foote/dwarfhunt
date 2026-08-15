"""dwarfhunt -- shared code for the brown dwarf / galaxy color-color work.

Nothing here imports species at module scope, so `import dwarfhunt` stays cheap
and side-effect free. Call `init()` to attach to the shared database.

The analysis modules -- planets, galaxies, gmm, plots -- do import species, so
they are deliberately NOT pulled in here. Import them explicitly:

    import dwarfhunt
    from dwarfhunt import planets, galaxies

    db = dwarfhunt.init()
"""

from . import paths
from .session import ensure_config, init

__all__ = ["init", "ensure_config", "paths"]
