import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


@pytest.fixture(scope="session")
def con():
    """Read-only connection to the built store. Integration tests skip without it."""
    import locus_db as db

    if not db.DB_PATH.exists():
        pytest.skip("locus.duckdb not built — run `uv run python ingest.py` first")
    # Read-only: tests never write, and a read-write handle would collide with
    # a running console holding the same file.
    connection = db.connect(read_only=True)
    yield connection
    connection.close()
