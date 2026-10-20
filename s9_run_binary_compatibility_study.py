"""이원 배합적합성 사례연구 — 1:1 (w/w) 기준.

s8_run_excipient_case_study.py는 선택한 부형제를 한꺼번에 엔진에 넣는
«전체 배합» 방식입니다. 그 방식은 부형제 사이에 중복 필터를 공유하기
때문에, 같은 불순물을 만드는 부형제가 여럿이면 먼저 처리된 하나만
남고 나머지는 «상호작용 없음»으로 보고됩니다.

이 스크립트는 부형제를 한 쌍씩 독립적으로 평가하고, 두 방식의 결과를
나란히 비교합니다. 비교표 자체가 논문에 쓸 수 있는 결과입니다.

2026-09-03 수정 — 비교 기준을 양쪽에 동일하게 적용 (C-6)
-------------------------------------------------------
이전 버전은 «전체 배합» 쪽을 엔진 원본 그대로(접힘 없음) 세고, «이원» 쪽은
`p17._collapse_positional_isomers` 로 **규칙당 최대 3개까지 접은 값**과
비교했습니다. 요약의 `(+N건)`은 서로 다른 두 집계 기준의 차이였고,
위치이성질체가 수십 건 나오는 폴리올 부형제에서는 오차가 수백 건 규모입니다.

이제 두 방식에 **같은 기준**을 적용합니다. `--basis` 로 고릅니다.
    collapsed (기본) : 양쪽 모두 p17 의 위치이성질체 접힘을 적용
    raw              : 양쪽 모두 접히지 않은 원본 건수
                       (이원 쪽 = len(pair.hits) + sum(pair.collapsed.values()))
어느 기준을 썼는지 리포트와 CSV에 반드시 기록됩니다.

2026-09-03 수정 — masking 판정을 양방향으로 (C-7)
-------------------------------------------------
이전 판정은 `in_matrix == 0` 만 봤기 때문에 실제로 일어난 재배분을 하나도
잡지 못했습니다. 이제 **건수 기준**(`전체배합 < 이원`)으로 판정하고,
반대 방향(전체배합이 더 많은 경우 = 과대·오귀속)도 따로 집계합니다.
p72 의 3단계 재니트로소화가 다른 부형제의 산물을 기질로 쓰기 때문에
전체배합 쪽에서 «아질산염 공급 부형제»로 몰아 귀속되는 건이 생깁니다
(p72 의 `substrate_origin_excipient` 필드로 확인합니다).

실행
    python s9_run_binary_compatibility_study.py
    python s9_run_binary_compatibility_study.py --all-excipients
    python s9_run_binary_compatibility_study.py --basis raw
    python s9_run_binary_compatibility_study.py Aspirin=CC(=O)Oc1ccccc1C(=O)O

출력
    data/binary_compatibility_report.txt
    data/binary_compatibility_pairs.csv
    data/binary_vs_matrix_comparison.csv

2026-08-25 작성.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

OUT_DIR = APP_DIR / "data"
OUT_TXT = OUT_DIR / "binary_compatibility_report.txt"
OUT_PAIRS = OUT_DIR / "binary_compatibility_pairs.csv"
OUT_COMPARE = OUT_DIR / "binary_vs_matrix_comparison.csv"

BAR = "=" * 78
_log: list[str] = []

# s8과 동일한 조합이라야 비교가 성립합니다
DEFAULT_EXCIPIENTS = [
    "Lactose_Monohydrate",
    "Microcrystalline_Cellulose",
    "Crospovidone",
    "Polyethylene_Glycol_4000",
    "Magnesium_Stearate",
]

DEFAULT_APIS: list[tuple[str, str]] = [
    ("Ranitidine", "CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1"),
    ("Atenolol", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"),
    ("Metformin", "CN(C)C(=N)NC(=N)N"),
]

RATIO = (1.0, 1.0)


def say(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


def main() -> int:
    raw_args = list(sys.argv[1:])
    use_all = "--all-excipients" in raw_args

    # 2026-09-03 추가(C-6): 비교 기준 선택 (--basis collapsed|raw, 기본 collapsed)
    basis = "collapsed"
    args: list[str] = []
    i = 0
    while i < len(raw_args):
        a = raw_args[i]
        if a == "--basis" and i + 1 < len(raw_args):
            basis = raw_args[i + 1].strip().lower()
            i += 2
            continue
        if a.startswith("--basis="):
            basis = a.split("=", 1)[1].strip().lower()
            i += 1
            continue
        if not a.startswith("--"):
            args.append(a)
        i += 1

    if basis not in ("collapsed", "raw"):
        print(f"[오류] --basis 는 collapsed 또는 raw 여야 합니다 (받은 값: {basis!r})")
        return 1

    BASIS_LABEL = {
        "collapsed": "위치이성질체 접힘 적용 (p17._collapse_positional_isomers, 규칙당 최대 3개) "
                     "— 양쪽 동일",
        "raw": "접힘 없는 원본 건수 — 양쪽 동일",
    }[basis]

    say(BAR)
    say("  PRIME-DS  Case Study — 이원 배합적합성 (binary drug–excipient, 1:1 w/w)")
    say(BAR)
    say()

    try:
        from p17_binary_compatibility import (
            _collapse_positional_isomers,
            evaluate_binary_compatibility,
        )
        from p72_excipient_interaction_engine import (
            EXCIPIENT_CHOICES,
            EXCIPIENT_LIBRARY,
            AdvancedExcipientMatrixEngine,
        )
    except Exception as exc:
        say(f"[오류] 엔진을 불러오지 못했습니다: {exc}")
        return 1

    # 2026-09-03 추가(C-6): 비교 기준을 리포트 첫머리에 못박는다.
    say("  [비교 기준]")
    say(f"    basis = {basis}   {BASIS_LABEL}")
    say("    이전 버전은 전체배합=원본 / 이원=접힘 으로 서로 다른 기준을 섞어")
    say("    비교했습니다. 그 «(+N건)» 은 방식 차이가 아니라 집계 기준 차이였습니다.")
    say()

    def count_matrix(payload: list[dict]) -> dict[str, int]:
        """전체 배합 산출물을 부형제별로 센다 — 이원 쪽과 같은 기준으로."""
        by_exc: dict[str, list[dict]] = {}
        for hit in payload:
            key = str(hit.get("excipient_name") or "")
            by_exc.setdefault(key, []).append(hit)
        if basis == "raw":
            return {k: len(v) for k, v in by_exc.items()}
        # collapsed — p17 과 동일한 함수를, p17 과 동일하게 «부형제 단위»로 적용
        out: dict[str, int] = {}
        for k, v in by_exc.items():
            kept, _ = _collapse_positional_isomers([dict(h) for h in v])
            out[k] = len(kept)
        return out

    def count_binary(pair) -> int:
        """이원 쪽 건수 — 같은 기준으로."""
        if basis == "raw":
            return len(pair.hits) + sum(int(x or 0) for x in (pair.collapsed or {}).values())
        return len(pair.hits)

    excipients = list(EXCIPIENT_CHOICES) if use_all else list(DEFAULT_EXCIPIENTS)

    apis = list(DEFAULT_APIS)
    if args:
        apis = []
        for i, raw in enumerate(args, start=1):
            if "=" in raw:
                name, smi = raw.split("=", 1)
                apis.append((name.strip(), smi.strip()))
            else:
                apis.append((f"API {i}", raw.strip()))

    say(f"배합비        : {RATIO[0]:g} : {RATIO[1]:g} (w/w, 원료 : 부형제)")
    say(f"부형제        : {len(EXCIPIENT_LIBRARY)}종 중 {len(excipients)}종")
    if not use_all:
        for key in excipients:
            entry = EXCIPIENT_LIBRARY.get(key, {})
            say(f"    · {entry.get('label_ko', key)}")
    say(f"원료의약품    : {len(apis)}종")
    say()

    pair_rows: list[dict] = []
    compare_rows: list[dict] = []

    for name, smi in apis:
        say(BAR)
        say(f"  {name}")
        say(BAR)

        # --- 전체 배합 (기존 방식) ---
        try:
            matrix_hits = AdvancedExcipientMatrixEngine().evaluate_formulation_risk(
                smi, excipients
            )
            matrix_payload = [
                h.to_dict() if hasattr(h, "to_dict") else dict(h) for h in (matrix_hits or [])
            ]
        except Exception as exc:
            say(f"  [오류] 전체 배합 평가 실패: {exc}")
            matrix_payload = []

        # 2026-09-03 수정(C-6): 이원 쪽과 동일한 기준으로 센다.
        matrix_by_exc = count_matrix(matrix_payload)
        matrix_total = sum(matrix_by_exc.values())

        # 2026-09-03 추가(C-7c): 전체배합의 오귀속을 보기 위한 «골격 제공 부형제» 집계.
        #   p72 의 3단계 재니트로소화는 다른 부형제의 산물을 기질로 쓰므로,
        #   excipient_name(= 아질산염 공급 부형제)만 보면 귀속이 한쪽으로 몰린다.
        origin_by_exc: dict[str, dict[str, int]] = {}
        for hit in matrix_payload:
            nitrite_exc = str(hit.get("excipient_name") or "")
            origin = str(hit.get("substrate_origin_excipient") or "")
            if origin and origin != nitrite_exc:
                origin_by_exc.setdefault(nitrite_exc, {})
                origin_by_exc[nitrite_exc][origin] = \
                    origin_by_exc[nitrite_exc].get(origin, 0) + 1

        # --- 이원 1:1 (신규) ---
        result = evaluate_binary_compatibility(smi, excipients, ratio=RATIO)
        if result.errors:
            for err in result.errors:
                say(f"  [오류] {err}")
            continue

        say(f"  SMILES        : {result.api_smiles}")
        say(f"  분자식 / MW   : {result.api_formula} · {result.api_mw:.2f} g/mol")
        say()
        binary_total = sum(count_binary(p) for p in result.pairs)
        say(f"  전체 배합 방식 : {matrix_total}건   (기준: {basis})")
        say(f"  이원 1:1 방식  : {binary_total}건   (기준: {basis})")
        say()

        say(f"    {'부형제':<32}{'점수':>7}{'구간':>10}{'M7':>5}{'건':>5}{'몰비':>8}"
            f"  {'전체배합':>8}{'Δ':>7}")
        say(f"    {'-' * 85}")
        masked: list[tuple[str, int]] = []       # 이원 > 전체배합 (전체배합이 놓친 것)
        over_attrib: list[tuple[str, int]] = []  # 전체배합 > 이원 (과대·오귀속 의심)
        for pair in result.pairs:
            in_matrix = matrix_by_exc.get(pair.excipient_key, 0)
            n_binary = count_binary(pair)
            delta = n_binary - in_matrix
            # 2026-09-03 수정(C-7a): `in_matrix == 0` 이 아니라 **건수 기준**으로 판정.
            #   0 만 보면 «7 → 4» 같은 실제 재배분을 하나도 잡지 못한다.
            if in_matrix < n_binary:
                masked.append((pair.excipient_key, delta))
            elif in_matrix > n_binary:
                over_attrib.append((pair.excipient_key, in_matrix - n_binary))
            mole = f"{pair.mole_ratio:.2f}" if pair.mole_ratio else "—"
            say(
                f"    {pair.excipient_key:<32}{pair.risk_score:>7.1f}{pair.risk_band:>10}"
                # 2026-09-03 수정: 확정 등급과 범위 추정("1-3")을 구분해 표시합니다.
                # 이전에는 범위 힌트에서 최소 자릿수를 뽑아 Class 1로 확정 표기했습니다.
                f"{str(getattr(pair, 'worst_m7_class_display', None) or pair.worst_m7_class or '-'):>5}"
                f"{n_binary:>5}{mole:>8}  {in_matrix:>8}{delta:>+7}"
            )

            for hit in pair.hits:
                pair_rows.append(
                    {
                        "API": name,
                        "API_SMILES": result.api_smiles,
                        "Excipient": pair.excipient_key,
                        "Excipient_MW": f"{pair.excipient_mw:.2f}" if pair.excipient_mw else "",
                        "Molar_ratio_at_1to1_ww": f"{pair.mole_ratio:.3f}" if pair.mole_ratio else "",
                        "Risk_score": f"{pair.risk_score:.1f}",
                        "Risk_band": pair.risk_band,
                        "Risk_class": hit.get("risk_class", ""),
                        "Rule": hit.get("rule_name", ""),
                        "ICH_M7_hint": hit.get("ich_m7_class_hint", ""),
                        "Alpha_hydrogen": hit.get("alpha_hydrogen_present", ""),
                        "Nitrite_exposure": hit.get("nitrite_exposure", ""),
                        "Impurity_SMILES": hit.get("smiles", ""),
                    }
                )

            # 2026-09-03 수정(C-7a): Delta·Matrix_only_hits 컬럼으로 **양방향** 보고
            compare_rows.append(
                {
                    "API": name,
                    "Excipient": pair.excipient_key,
                    "Comparison_basis": basis,
                    "Matrix_mode_hits": in_matrix,
                    "Binary_mode_hits": n_binary,
                    "Delta_binary_minus_matrix": delta,
                    "Masked_in_matrix_mode": "YES" if in_matrix < n_binary else "",
                    "Matrix_only_hits": max(0, in_matrix - n_binary),
                    "Matrix_substrate_origin_mismatch": sum(
                        origin_by_exc.get(pair.excipient_key, {}).values()
                    ),
                    "Risk_score": f"{pair.risk_score:.1f}",
                    "Risk_band": pair.risk_band,
                }
            )

        say()
        # 2026-09-03 수정(C-7a): 두 방향을 모두 보고한다.
        if masked:
            say(f"  ▽ 전체 배합 방식이 **적게** 보고한 부형제 {len(masked)}종 (이원 > 전체배합)")
            for key, d in masked:
                entry = EXCIPIENT_LIBRARY.get(key, {})
                say(f"        · {entry.get('label_ko', key)}  ({d:+}건)")
            say("      같은 불순물을 만드는 다른 부형제가 먼저 처리되어")
            say("      건수가 줄거나 «상호작용 없음»으로 보고되고 있었습니다.")
        if over_attrib:
            say(f"  △ 전체 배합 방식이 **많게** 보고한 부형제 {len(over_attrib)}종 "
                "(전체배합 > 이원 · 과대/오귀속 의심)")
            for key, d in over_attrib:
                entry = EXCIPIENT_LIBRARY.get(key, {})
                mism = sum(origin_by_exc.get(key, {}).values())
                extra = ""
                if mism:
                    top = sorted(origin_by_exc.get(key, {}).items(), key=lambda kv: -kv[1])
                    extra = ("   ← 골격 제공 부형제가 다른 건 "
                             f"{mism}건 ({', '.join(f'{k}:{v}' for k, v in top[:3])})")
                say(f"        · {entry.get('label_ko', key)}  (+{d}건){extra}")
            say("      p72 의 3단계 재니트로소화는 **다른 부형제의 산물**을 기질로 쓰는데,")
            say("      결과는 아질산염을 공급한 부형제에 귀속됩니다. 전체배합에서만")
            say("      성립하는 경로이므로 이원 평가에는 나타나지 않습니다.")
            say("      → 이 차이는 «위음성»이 아니라 «귀속 위치 문제»입니다.")
        if not masked and not over_attrib:
            say("  두 방식의 부형제별 건수가 모두 일치합니다.")
        say()
        clean = [p.excipient_key for p in result.pairs if not p.hits]
        if clean:
            say(f"  이상 없음 : {', '.join(clean)}")
        say()

    # ------------------------------------------------------------------
    say(BAR)
    say("  요약")
    say(BAR)
    total_matrix = sum(int(r["Matrix_mode_hits"]) for r in compare_rows)
    total_binary = sum(int(r["Binary_mode_hits"]) for r in compare_rows)
    total_masked = sum(1 for r in compare_rows if r["Masked_in_matrix_mode"])
    total_over = sum(1 for r in compare_rows if int(r["Matrix_only_hits"]) > 0)
    sum_masked = sum(max(0, int(r["Delta_binary_minus_matrix"])) for r in compare_rows)
    sum_over = sum(int(r["Matrix_only_hits"]) for r in compare_rows)
    sum_mismatch = sum(int(r["Matrix_substrate_origin_mismatch"]) for r in compare_rows)
    say(f"  비교 기준        : {basis}  —  {BASIS_LABEL}")
    say(f"  평가 쌍          : {len(compare_rows)}")
    say(f"  전체 배합 방식   : {total_matrix}건")
    say(f"  이원 1:1 방식    : {total_binary}건  ({total_binary - total_matrix:+}건)")
    say(f"  이원이 더 많음   : {total_masked}쌍 (합계 {sum_masked:+}건)  ← 전체배합의 과소보고")
    say(f"  전체배합이 더 많음: {total_over}쌍 (합계 +{sum_over}건)  ← 과대/오귀속 의심")
    say(f"  골격 제공 부형제가 다른 건 : {sum_mismatch}건 "
        "(p72 3단계 재니트로소화 — 아질산염 공급자에 귀속됨)")
    say()

    # 2026-09-03 수정(C-7b): 결론을 관측치로 가드한다.
    #   이전에는 masking 관측이 0건이어도 «matrix 방식은 위음성을 만든다» 는
    #   단정형 결론을 무조건 인쇄했다.
    if total_masked > 0:
        say("  [관측] 전체 배합 방식은 부형제 사이에 중복 필터를 공유하므로, 같은")
        say("  불순물을 만드는 부형제가 여럿이면 일부만 보고됩니다. 선택 순서에")
        say(f"  따라 결과가 달라지며, 이 실행에서는 {total_masked}쌍 · 합계 {sum_masked}건이")
        say("  전체배합에서 과소보고되었습니다 (위양성이 아니라 위음성).")
    else:
        say("  [관측] 이 조합에서는 masking(전체배합의 과소보고)이 관측되지 않았습니다.")
    if total_over > 0:
        say()
        say(f"  [관측] 반대 방향도 있습니다 — {total_over}쌍 · 합계 +{sum_over}건에서")
        say("  전체 배합 방식이 **더 많이** 보고했습니다. 이는 과소보고가 아니라")
        say("  과대·오귀속입니다. p72 의 3단계 재니트로소화가 다른 부형제의 산물을")
        say("  기질로 쓰는데 결과는 아질산염 공급 부형제에 귀속되기 때문이며,")
        say("  전체배합에서만 성립하는 경로라 이원 평가에는 나타나지 않습니다.")
        say("  두 방식의 차이를 «위음성»으로만 해석하면 방향이 반대가 됩니다.")
    say()
    say("  ※ 위 수치는 모두 basis="
        f"{basis} 기준입니다. 다른 기준으로 보려면 --basis 로 바꿔 실행하십시오.")
    say()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        OUT_TXT.write_text("\n".join(_log), encoding="utf-8")
        if pair_rows:
            with OUT_PAIRS.open("w", newline="", encoding="utf-8-sig") as fh:
                writer = csv.DictWriter(fh, fieldnames=list(pair_rows[0].keys()))
                writer.writeheader()
                writer.writerows(pair_rows)
        if compare_rows:
            with OUT_COMPARE.open("w", newline="", encoding="utf-8-sig") as fh:
                writer = csv.DictWriter(fh, fieldnames=list(compare_rows[0].keys()))
                writer.writeheader()
                writer.writerows(compare_rows)
        print(f"요약 기록 : data/{OUT_TXT.name}")
        print(f"쌍별 결과 : data/{OUT_PAIRS.name}")
        print(f"방식 비교 : data/{OUT_COMPARE.name}")
    except Exception as exc:
        print(f"(저장 실패: {exc})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
