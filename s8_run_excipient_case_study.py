"""
논문 Case Study용 — 부형제-API 상호작용 예측 (완제의약품/DP).

원료의약품과 부형제를 조합해 생성 가능한 불순물을 예측하고,
ICH M7 4층 앙상블로 판정합니다. 니트로사민 생성물에는
run_case_study.py 와 동일한 검증 규칙(순전하 0 · N-N=O 구조)을 적용합니다.

사용법
------
    python run_excipient_case_study.py
    python run_excipient_case_study.py 이름=SMILES

결과물
------
    data/excipient_case_study_results.csv
    data/excipient_case_study_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUT_CSV = DATA_DIR / "excipient_case_study_results.csv"
OUT_TXT = DATA_DIR / "excipient_case_study_report.txt"

STRICT_NITROSAMINE_SMARTS = "[#7]-[#7X2;+0]=[OX1;+0]"

# 아질산염 잔류 프로파일이 등록된 부형제 위주로 구성
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

BAR = "=" * 78
_log: list[str] = []


def say(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


def main() -> int:
    say(BAR)
    say("  PRIME-DS  Case Study — 부형제-API 상호작용 (완제의약품)")
    say(BAR)
    say()

    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors

    RDLogger.DisableLog("rdApp.*")

    from p72_excipient_interaction_engine import (
        EXCIPIENT_LIBRARY,
        AdvancedExcipientMatrixEngine,
    )
    from p65_m7_ensemble_matrix import run_m7_ensemble_matrix

    try:
        from p60_m7_empirical_db import get_default_db_connection
        conn = get_default_db_connection()
    except Exception as exc:  # noqa: BLE001
        say(f"(실측 DB 연결 실패: {exc} — Layer A 없이 진행)")
        conn = None

    import p63_deepchem_predictor as dp
    say(f"Layer C 모델 : {dp.model_status()['model_label']}")
    say(f"부형제 라이브러리 : {len(EXCIPIENT_LIBRARY)}종 중 {len(DEFAULT_EXCIPIENTS)}종 선택")
    for k in DEFAULT_EXCIPIENTS:
        entry = EXCIPIENT_LIBRARY.get(k, {})
        say(f"    · {entry.get('label_ko', k)}")
    say()

    args = sys.argv[1:]
    apis = DEFAULT_APIS
    if args:
        apis = []
        for i, a in enumerate(args, start=1):
            if "=" in a:
                n, s = a.split("=", 1)
                apis.append((n.strip(), s.strip()))
            else:
                apis.append((f"API {i}", a.strip()))

    engine = AdvancedExcipientMatrixEngine()
    strict_q = Chem.MolFromSmarts(STRICT_NITROSAMINE_SMARTS)
    rows: list[dict] = []

    for name, api_smi in apis:
        say(BAR)
        say(f"  {name}")
        say(BAR)
        mol = Chem.MolFromSmiles(api_smi)
        if mol is None:
            say(f"  [오류] SMILES 파싱 실패: {api_smi}")
            continue
        say(f"  SMILES : {Chem.MolToSmiles(mol)}")

        try:
            hits = engine.evaluate_formulation_risk(api_smi, DEFAULT_EXCIPIENTS)
        except Exception as exc:  # noqa: BLE001
            say(f"  [오류] 예측 실패: {exc}")
            continue

        say(f"  원시 히트 : {len(hits)}건")

        seen: dict[str, int] = {}
        kept: list[dict] = []
        n_charge = n_bad_nitroso = n_dup = 0

        for h in hits:
            d = h.to_dict() if hasattr(h, "to_dict") else dict(h)
            smi = str(d.get("smiles") or "")
            m2 = Chem.MolFromSmiles(smi) if smi else None
            if m2 is None:
                continue

            rule = str(d.get("rule_name") or "")
            is_nitroso_rule = "ndsri" in rule.lower() or "nitros" in rule.lower()

            # 전하 불균형은 어떤 경로든 표기 오류로 간주
            if Chem.GetFormalCharge(m2) != 0:
                n_charge += 1
                continue
            # 니트로사민 경로 생성물은 N-N=O 구조를 반드시 가져야 함
            if is_nitroso_rule and not m2.HasSubstructMatch(strict_q):
                n_bad_nitroso += 1
                continue

            key = Chem.MolToInchiKey(m2)
            if key in seen:
                n_dup += 1
                continue
            seen[key] = 1

            canon = Chem.MolToSmiles(m2)
            try:
                res = run_m7_ensemble_matrix(canon, db_connection=conn)
                cls = res.class_label
                prob = res.layer_c.p_mutagenic
                alerts = "; ".join(getattr(res.layer_b, "combined_alerts", []) or [])
            except Exception:  # noqa: BLE001
                cls, prob, alerts = "판정 실패", None, ""

            # 2026-09-15 추가 — 부형제 귀속 두 축을 함께 내보낸다.
            #
            #   excipient (= excipient_name) 는 하위 호환 열이며 3단계
            #   (산물 재니트로소화)에서는 **아질산염을 공급한 부형제**만
            #   가리킨다. 반응 골격을 준 부형제는 다른 것일 수 있다.
            #   p72 가 2026-09-03 에 두 필드를 추가했는데 s8 이 내보내지
            #   않아, CSV 만 보면 귀속을 확인할 수 없었다.
            #
            #   ⚠ 어느 열로 집계하든 «부형제별 순위»는 s8 로 내지 말 것.
            #      evaluate_formulation_risk 는 seen_smiles 를 부형제
            #      사이에 공유하므로 같은 산물을 내는 부형제가 여럿이면
            #      먼저 처리된 하나만 남는다. 부형제 단위 비교는 부형제마다
            #      엔진을 새로 만드는 s9(이원 1:1)를 쓴다.
            kept.append({
                "excipient": d.get("excipient_name", ""),
                "nitrite_source": d.get("nitrite_source_excipient", ""),
                "substrate_origin": d.get("substrate_origin_excipient", ""),
                "rule": rule,
                "smiles": canon,
                "mw": round(Descriptors.MolWt(m2), 2),
                "ich_m7": cls,
                "prob": prob,
                "alerts": alerts,
                "alpha_h": d.get("alpha_hydrogen_present"),
                "steric": d.get("steric_hindrance"),
                "nitrite": d.get("nitrite_exposure"),
                "confidence": d.get("confidence", ""),
            })

        if n_charge or n_bad_nitroso or n_dup:
            say(f"    제외 — 전하 표기 오류 {n_charge} · "
                f"N-니트로사민 아님 {n_bad_nitroso} · 중복 {n_dup}")
        say(f"  검증된 불순물 : {len(kept)}건")
        say()

        if kept:
            say(f"    {'부형제':<22} {'기전':<24} {'ICH M7':<9} {'P':>6}  MW")
            say("    " + "-" * 72)
            for k in kept:
                p = f"{k['prob']:.3f}" if isinstance(k["prob"], float) else "  n/a"
                say(f"    {str(k['excipient'])[:20]:<22} {k['rule'][:22]:<24} "
                    f"{str(k['ich_m7'])[:8]:<9} {p:>6}  {k['mw']:.1f}")
                say(f"        {k['smiles']}")
                if k.get("nitrite"):
                    say(f"        아질산염 {k['nitrite']} · alpha-H {k['alpha_h']} · 입체장애 {k['steric']}")
            say()

        for k in kept:
            k2 = dict(k)
            k2["api_name"] = name
            k2["api_smiles"] = api_smi
            rows.append(k2)

    # 요약
    say(BAR)
    say("  요약")
    say(BAR)
    say(f"  검토 원료의약품 : {len(apis)}종 × 부형제 {len(DEFAULT_EXCIPIENTS)}종")
    say(f"  검증된 불순물   : {len(rows)}건")
    if rows:
        from collections import Counter
        say()
        say("  기전별")
        for r, n in Counter(x["rule"] for x in rows).most_common():
            say(f"      {r}: {n}건")
        say()
        say("  ICH M7 등급별")
        for c, n in Counter(str(x["ich_m7"]) for x in rows).most_common():
            say(f"      {c}: {n}건")
        say()
        say("  부형제별 — ⚠ 논문 인용 금지")
        say("    evaluate_formulation_risk 는 중복 제거 집합을 부형제 사이에")
        say("    공유한다. 같은 산물을 내는 부형제가 여럿이면 먼저 처리된")
        say("    하나만 남으므로 아래 분포는 DEFAULT_EXCIPIENTS 의 순서에")
        say("    의존하는 인공물이다. 부형제 단위 비교는 s9(이원 1:1)를 쓸 것.")
        say("    순서에 불변인 것은 총 건수 · 기전 분포 · 등급 분포뿐이다.")
        say()
        say("    아질산염 공급 부형제 기준 (excipient_name)")
        for e, n in Counter(str(x["excipient"]) for x in rows).most_common():
            say(f"      {e}: {n}건")

        # 3단계에서는 골격을 준 부형제가 아질산염 공급자와 다를 수 있다.
        subs = Counter(str(x.get("substrate_origin") or "") for x in rows)
        subs.pop("", None)
        if subs:
            say()
            say("    반응 골격 제공 부형제 기준 (substrate_origin_excipient)")
            for e, n in subs.most_common():
                say(f"      {e}: {n}건")
            say()
            say("    두 분포가 다르면 3단계 산물의 귀속이 두 부형제에")
            say("    걸쳐 있다는 뜻이다.")
    say()

    try:
        import pandas as pd
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
        OUT_TXT.write_text("\n".join(_log), encoding="utf-8")
        print(f"상세 결과: data/{OUT_CSV.name}")
        print(f"요약 기록: data/{OUT_TXT.name}")
    except Exception as exc:  # noqa: BLE001
        print(f"(저장 실패: {exc})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
