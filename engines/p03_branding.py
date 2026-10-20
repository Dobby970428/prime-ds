"""PRIME-DS Process Mode — UI 브랜딩 · 터미널 로그 접두사."""

from __future__ import annotations

import logging

LOG_PREFIX = "[PRIME-DS-PROCESS]"
UI_ENGINE_TITLE = "PRIME-DS 공정 조건 예측 엔진"
UI_PROCESS_MODE_TITLE = "PRIME-DS Process Mode"
UI_RESULTS_TITLE = f"{UI_ENGINE_TITLE} — 예측 결과"

_PROCESS_LOGGER = logging.getLogger("prime_ds.process")


def get_process_logger() -> logging.Logger:
    if not _PROCESS_LOGGER.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(message)s"))
        _PROCESS_LOGGER.addHandler(handler)
        _PROCESS_LOGGER.setLevel(logging.INFO)
    return _PROCESS_LOGGER


def log_process_info(message: str) -> None:
    get_process_logger().info(f"{LOG_PREFIX} {message}")


def log_process_warning(message: str) -> None:
    get_process_logger().warning(f"{LOG_PREFIX} {message}")


def log_process_error(message: str) -> None:
    get_process_logger().error(f"{LOG_PREFIX} {message}")


def log_process_exception(message: str, exc: BaseException | None = None) -> None:
    get_process_logger().exception(f"{LOG_PREFIX} {message}", exc_info=exc)
