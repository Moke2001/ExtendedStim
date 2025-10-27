"""ExtendStim package initializer.
Exports the main subpackages and registers a lowercase alias so users can import
either `ExtendStim` or `extendstim` on case-insensitive filesystems (Windows).
"""

# Re-export subpackages for convenience
from . import Circuit, Code, Math, Physics, Platform  # noqa: F401

__all__ = ["Circuit", "Code", "Math", "Physics", "Platform"]

# Package version (updated by releases)
__version__ = "0.1.0"

# Ensure lowercase alias exists in sys.modules so `import extendstim` works
import sys
# Register both canonical and lowercase names to the same module object
mod = sys.modules.get(__name__)
if mod is not None:
    # Do not overwrite existing mappings if they point somewhere else
    if 'extendstim' not in sys.modules:
        sys.modules['extendstim'] = mod
    if 'ExtendStim' not in sys.modules:
        sys.modules['ExtendStim'] = mod