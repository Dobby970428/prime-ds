"""Build or rebuild the ICH M7 empirical SQLite index from seed CSV."""

from __future__ import annotations

from p60_m7_empirical_db import build_m7_empirical_db, empirical_db_stats

if __name__ == "__main__":
    count = build_m7_empirical_db(force=True)
    stats = empirical_db_stats()
    print(f"Built {count} records -> {stats}")
