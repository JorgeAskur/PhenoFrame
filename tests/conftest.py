"""Shared pytest fixtures and helpers for the phenosuite test suite."""

from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SAMPLE_PLANT_XML = REPO_ROOT / "plants" / "plant_0.xml"


def _wrapper_available() -> bool:
    """Return True when the C++ shared library can actually be loaded."""
    try:
        from phenosuite import Maize  # noqa: WPS433 (intentional local import)

        with Maize():
            return True
    except (FileNotFoundError, OSError, RuntimeError):
        return False


# Evaluated once at import time; used by `requires_wrapper` to skip cleanly
# on platforms without a compiled shared library (e.g. CI before build step).
WRAPPER_AVAILABLE = _wrapper_available()

requires_wrapper = pytest.mark.skipif(
    not WRAPPER_AVAILABLE,
    reason="Compiled maize_c_api shared library not available on this platform",
)


@pytest.fixture(scope="session")
def sample_plant_xml() -> Path:
    """Return the path to the bundled sample descriptor."""
    if not SAMPLE_PLANT_XML.exists():
        pytest.skip(f"Sample descriptor missing: {SAMPLE_PLANT_XML}")
    return SAMPLE_PLANT_XML
