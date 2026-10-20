"""분해산물 엔진·UI session 캐시 버전 (순환 import 방지용 단독 상수)."""

# 상수 하나

# 엔진 로직 변경 시 반드시 증가 — 구 Mock(DP-) 캐시 자동 무효화
ENGINE_CACHE_VERSION = "v16-toxic-degradation-tree"
