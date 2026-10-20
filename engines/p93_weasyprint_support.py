"""PRIME-DS — WeasyPrint Windows GTK/Pango DLL 경로 설정."""

from __future__ import annotations

import os
import sys
from pathlib import Path

_DEFAULT_DLL_DIRS: tuple[Path, ...] = (
    Path(r"C:\msys64\mingw64\bin"),
    Path("/mingw64/bin"),
)


def configure_weasyprint_dll_paths() -> str | None:
    """WEASYPRINT_DLL_DIRECTORIES 및 Windows DLL search path 설정."""
    existing = (os.environ.get("WEASYPRINT_DLL_DIRECTORIES") or "").strip()
    if existing and Path(existing).is_dir():
        _register_dll_dir(existing)
        return existing

    for candidate in _DEFAULT_DLL_DIRS:
        if candidate.is_dir():
            path_str = str(candidate.resolve())
            os.environ["WEASYPRINT_DLL_DIRECTORIES"] = path_str
            _register_dll_dir(path_str)
            return path_str
    return None


def _register_dll_dir(path_str: str) -> None:
    if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
        try:
            os.add_dll_directory(path_str)
        except OSError:
            pass


def weasyprint_ready() -> bool:
    """WeasyPrint import 가능 여부."""
    configure_weasyprint_dll_paths()
    try:
        from weasyprint import HTML  # noqa: F401

        return True
    except Exception:
        return False
