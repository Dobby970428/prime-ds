"""PRIME-DS Process Mode — 연산 세이프가드 · 타임아웃 · 제품 수 상한."""

from __future__ import annotations

import time
from typing import Any

from engines.p03_branding import log_process_warning

MAX_PROCESS_PRODUCTS: int = 500
PROCESS_REQUEST_TIMEOUT_SEC: float = 15.0
INVALID_SMILES_BYPASS_MSG = "Invalid SMILES or Valency Error bypassed."


def log_invalid_smiles_bypass(*, context: str = "") -> None:
    """파싱·원자가 오류 시 표준 bypass 로그."""
    if context:
        log_process_warning(f"{INVALID_SMILES_BYPASS_MSG} context={context[:120]}")
    else:
        log_process_warning(INVALID_SMILES_BYPASS_MSG)


class ProcessRunContext:
    """단일 예측 요청 내 제품 수 상한 · 경과 시간 제한."""

    __slots__ = ("deadline", "max_products", "timeout_sec", "product_count", "timed_out", "capped")

    def __init__(
        self,
        *,
        max_products: int = MAX_PROCESS_PRODUCTS,
        timeout_sec: float = PROCESS_REQUEST_TIMEOUT_SEC,
        autostart: bool = True,
    ) -> None:
        self.max_products = max_products
        self.timeout_sec = timeout_sec
        self.product_count = 0
        self.timed_out = False
        self.capped = False
        self.deadline = 0.0
        if autostart:
            self.restart_timer()

    def restart_timer(self) -> None:
        """연산 구간 타이머 재시작 (캐시 워밍 등 제외)."""
        self.deadline = time.monotonic() + self.timeout_sec
        self.timed_out = False

    def check_continue(self) -> bool:
        if time.monotonic() >= self.deadline:
            if not self.timed_out:
                self.timed_out = True
                log_process_warning(
                    f"Request timeout ({PROCESS_REQUEST_TIMEOUT_SEC:.0f}s); "
                    "computation truncated."
                )
            return False
        if self.product_count >= self.max_products:
            if not self.capped:
                self.capped = True
                log_process_warning(
                    f"Product cap reached ({self.max_products}); "
                    "remaining candidates skipped."
                )
            return False
        return True

    def register_product(self) -> bool:
        if not self.check_continue():
            return False
        self.product_count += 1
        return True

    def truncate_records(self, records: list[dict[str, Any]]) -> list[dict[str, Any]]:
        if len(records) <= self.max_products:
            return records
        self.capped = True
        log_process_warning(
            f"Product cap reached ({self.max_products}); output truncated."
        )
        return records[: self.max_products]
