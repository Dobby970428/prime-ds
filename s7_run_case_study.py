"""
논문 Case Study용 — 니트로사민(NDSRI) 생성 예측 + ICH M7 판정.

지정한 원료의약품에 아질산염 조건을 적용해 생성 가능한 니트로사민 구조를
예측하고, 각각을 4층 앙상블(Layer A/B/C/D)로 판정해 표로 출력합니다.

사용법
------
    python run_case_study.py                       # 기본 4종 (메트포민·발사르탄·라니티딘·니자티딘)
    python run_case_study.py "CN(C)C(=N)NC(=N)N"   # SMILES 직접 지정
    python run_case_study.py 이름=SMILES 이름2=SMILES2

결과물
------
    data/case_study_results.csv    화합물별 상세
    data/case_study_report.txt     화면 출력 그대로
"""

from __future__ import annotations

import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUT_CSV = DATA_DIR / "case_study_results.csv"
OUT_TXT = DATA_DIR / "case_study_report.txt"

BAR = "=" * 76
_log: list[str] = []

# 생성물이 진짜 N-니트로사민인지 검증하는 패턴.
#   [#7]      아민 질소
#   [#7X2;+0] 니트로소 질소 — 결합 2개, 전하 0 (전하가 붙어 있으면 표기 오류)
#   [OX1;+0]  말단 산소, 전하 0
# C-니트로소(탄소에 붙은 N=O)는 이 패턴에 걸리지 않는다.
STRICT_NITROSAMINE_SMARTS = "[#7]-[#7X2;+0]=[OX1;+0]"


def say(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


# 기본 사례 — 니트로사민 이슈로 실제 회수가 있었던 품목 위주
DEFAULT_CASES: list[tuple[str, str]] = [
    ("Metformin", "CN(C)C(=N)NC(=N)N"),
    ("Valsartan", "CCCCC(=O)N(Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1)C(C(C)C)C(=O)O"),
    ("Ranitidine", "CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1"),
    ("Nizatidine", "CNC(=C[N+](=O)[O-])NCCSCc1csc(CN(C)C)n1"),
]


def parse_args() -> list[tuple[str, str]]:
    args = sys.argv[1:]
    if not args:
        return DEFAULT_CASES
    cases: list[tuple[str, str]] = []
    for i, a in enumerate(args, start=1):
        if "=" in a:
            name, smi = a.split("=", 1)
            cases.append((name.strip(), smi.strip()))
        else:
            cases.append((f"Compound {i}", a.strip()))
    return cases


def main() -> int:
    say(BAR)
    say("  PRIME-DS  Case Study — 니트로사민(NDSRI) 생성 예측 및 ICH M7 판정")
    say(BAR)
    say()

    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Descriptors, rdMolDescriptors
    except ImportError as exc:
        say(f"[오류] RDKit 없음: {exc}")
        return 1
    RDLogger.DisableLog("rdApp.*")

    from p30_reaction_engine import (
        NITRITE_SMILES,
        detect_amine_sites,
        predict_nitrosation_with_reagent,
    )

    # 판정 엔진
    try:
        from p60_m7_empirical_db import get_default_db_connection
        conn = get_default_db_connection()
    except Exception as exc:  # noqa: BLE001
        say(f"(실측 DB 연결 실패: {exc} — Layer A 없이 진행)")
        conn = None
    from p65_m7_ensemble_matrix import run_m7_ensemble_matrix

    import p63_deepchem_predictor as dp
    st = dp.model_status()
    say(f"Layer C 모델 : {st['model_label']}")
    say(f"아질산염 시약 : {NITRITE_SMILES}")
    say()

    cases = parse_args()
    rows: list[dict] = []

    for name, api_smi in cases:
        say(BAR)
        say(f"  {name}")
        say(BAR)

        mol = Chem.MolFromSmiles(api_smi)
        if mol is None:
            say(f"  [오류] SMILES를 읽을 수 없습니다: {api_smi}")
            say()
            continue

        say(f"  SMILES : {Chem.MolToSmiles(mol)}")
        say(f"  분자량  : {Descriptors.MolWt(mol):.2f}")

        info = detect_amine_sites(api_smi)
        say(f"  니트로소화 가능 아민 : "
            f"{'있음' if info.has_nitrosatable_amine else '없음'}"
            f"  (2급 {getattr(info, 'secondary_count', '?')} · "
            f"3급 {getattr(info, 'tertiary_count', '?')})")

        if not info.has_nitrosatable_amine:
            say("  → 니트로소화 경로 없음. 니트로사민 생성 예측 없음.")
            say()
            rows.append({"api_name": name, "api_smiles": api_smi,
                         "nitrosatable": False, "product_smiles": "",
                         "final_class": "", "layer_b_alerts": "",
                         "layer_c_prob": "", "note": "no nitrosatable amine"})
            continue

        _, product_mols, _ = predict_nitrosation_with_reagent(api_smi, NITRITE_SMILES)

        # --- 검증 + 중복 제거 -------------------------------------------
        # 반응 규칙이 내놓는 원시 생성물에는 (a) 전하 표기가 깨진 구조와
        # (b) 탄소에 니트로소가 붙은 C-니트로소 화합물이 섞여 있다.
        # 둘 다 NDSRI가 아니므로 걸러내고, InChIKey로 중복을 제거한다.
        strict_q = Chem.MolFromSmarts(STRICT_NITROSAMINE_SMARTS)
        seen_keys: dict[str, str] = {}
        products: list[str] = []
        n_raw = 0
        n_charge = 0
        n_not_nitrosamine = 0
        n_dup = 0

        for pm in product_mols:
            try:
                s = Chem.MolToSmiles(pm)
                m2 = Chem.MolFromSmiles(s)
            except Exception:  # noqa: BLE001
                continue
            if not s or m2 is None:
                continue
            n_raw += 1

            if Chem.GetFormalCharge(m2) != 0:
                n_charge += 1          # 전하 불균형 → 표기 오류
                continue
            if not m2.HasSubstructMatch(strict_q):
                n_not_nitrosamine += 1  # N-N=O 없음 → 니트로사민 아님
                continue

            key = Chem.MolToInchiKey(m2)
            if key in seen_keys:
                n_dup += 1
                continue
            seen_keys[key] = s
            products.append(s)

        say(f"  원시 생성물      : {n_raw}종")
        if n_charge or n_not_nitrosamine or n_dup:
            say(f"    제외 — 전하 표기 오류 {n_charge} · "
                f"N-니트로사민 아님 {n_not_nitrosamine} · 중복 {n_dup}")
        say(f"  검증된 니트로사민 : {len(products)}종")
        say()

        if not products:
            rows.append({"api_name": name, "api_smiles": api_smi,
                         "nitrosatable": True, "product_smiles": "",
                         "final_class": "", "layer_b_alerts": "",
                         "layer_c_prob": "", "note": "no product generated"})
            continue

        say(f"    {'#':<3} {'ICH M7':<10} {'P(변이원성)':>11}  {'MW':>8}  구조")
        say("    " + "-" * 70)

        for idx, psmi in enumerate(products, start=1):
            try:
                res = run_m7_ensemble_matrix(psmi, db_connection=conn)
            except Exception as exc:  # noqa: BLE001
                say(f"    {idx:<3} 판정 실패: {exc}")
                continue

            pmol = Chem.MolFromSmiles(psmi)
            mw = Descriptors.MolWt(pmol) if pmol else 0.0
            formula = rdMolDescriptors.CalcMolFormula(pmol) if pmol else ""
            prob = res.layer_c.p_mutagenic
            prob_s = f"{prob:.3f}" if prob is not None else "  n/a"
            alerts = "; ".join(getattr(res.layer_b, "combined_alerts", []) or [])

            say(f"    {idx:<3} {res.class_label:<10} {prob_s:>11}  {mw:8.2f}  {psmi}")
            if alerts:
                say(f"        경보: {alerts}")

            rows.append({
                "api_name": name,
                "api_smiles": api_smi,
                "nitrosatable": True,
                "product_no": idx,
                "product_smiles": psmi,
                "product_formula": formula,
                "product_mw": round(mw, 2),
                "final_class": res.class_label,
                "rule_class": res.rule_class,
                "layer_a_db_hit": bool(res.layer_a.is_in_db),
                "layer_b_alert": bool(res.layer_b.alert_present),
                "layer_b_alerts": alerts,
                "layer_c_prob": prob,
                "ai_inferred": bool(res.ai_inferred_liability),
                "note": "",
            })
        say()

    # 요약
    say(BAR)
    say("  요약")
    say(BAR)
    made = [r for r in rows if r.get("product_smiles")]
    say(f"  검토한 원료의약품 : {len(cases)}종")
    say(f"  생성된 니트로사민 : {len(made)}종")
    if made:
        from collections import Counter
        dist = Counter(r["final_class"] for r in made)
        for cls, n in sorted(dist.items()):
            say(f"      {cls}: {n}종")
        alerted = sum(1 for r in made if r.get("layer_b_alert"))
        say(f"  구조 경보 양성   : {alerted}/{len(made)}종")
        probs = [r["layer_c_prob"] for r in made if isinstance(r.get("layer_c_prob"), float)]
        if probs:
            say(f"  Layer C 확률 범위 : {min(probs):.3f} ~ {max(probs):.3f}")
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
