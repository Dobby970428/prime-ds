"""
논문 Case Study (DS) — 원료의약품 공정 유래 불순물 생성 + ICH M7 판정.

로사르탄 상용 공정의 두 단계를 출발물질만으로 재현한다.

  1단계  Suzuki-Miyaura 커플링 (OTBN 생성)
         p36_suzuki_impurity_rules 의 규칙 5종
  2단계  테트라졸 형성 단계의 아자이드 친핵 치환
         p31_process_reaction_rules 의 Nucleophilic_Substitution_with_Azide

생성된 구조마다 물성(MW·LogP·끓는점 추정) → 4계층 ICH M7 → 적용범위(d_max)
전 과정을 돌리고, 문헌·규제 문서에 실제로 보고된 공정 불순물과 대조한다.

결과물
------
    data/ds_process_case_study_results.csv
    data/ds_process_case_study_report.txt
"""

from __future__ import annotations

import csv
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
OUT_CSV = DATA_DIR / "ds_process_case_study_results.csv"
OUT_TXT = DATA_DIR / "ds_process_case_study_report.txt"

BAR = "=" * 78
_log: list[str] = []


def say(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


# ---------------------------------------------------------------------------
# 공정 정의 — 로사르탄
# ---------------------------------------------------------------------------
SUZUKI_HALIDE = ("2-bromobenzonitrile", "N#Cc1ccccc1Br")
SUZUKI_BORONIC = ("4-methylphenylboronic acid", "Cc1ccc(B(O)O)cc1")

AZIDE = ("azide (NaN3 / HN3)", "[N-]=[N+]=N")
AZIDE_SUBSTRATES = [
    ("losartan chloromethyl intermediate",
     "CCCCc1nc(Cl)c(CCl)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1"),
    ("2-cyano-4'-(bromomethyl)biphenyl tetrazole (BBB-tetrazole)",
     "BrCc1ccc(-c2ccccc2-c2nn[nH]n2)cc1"),
]

LOSARTAN = "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1"

# 문헌·규제 문서에 실제로 보고된 공정 불순물 (구조로 대조한다)
#   key = 정규 SMILES
REPORTED: dict[str, dict[str, str]] = {
    "N#Cc1ccccc1-c1ccccc1C#N": {
        "name": "2,2'-dicyanobiphenyl",
        "cas": "4341-02-0",
        "source": "US6153810A — o-halobenzonitrile 호모커플링 부산물, 실시예에서 0.45% · 1.13% 로 정량",
    },
    "Cc1ccc(-c2ccc(C)cc2)cc1": {
        "name": "4,4'-dimethylbiphenyl",
        "cas": "613-33-2",
        "source": "US6392080 — p-tolyl 유기금속 자체 축합 부산물, 종래법 8–12% → 개선법 0.5–1%",
    },
    "CCCCc1nc(Cl)c(CN=[N+]=[N-])n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1": {
        "name": "losartan azido impurity",
        "cas": "727718-93-6",
        "source": "EDQM 2021 — Ames 양성, ICH M7 TTC 이하 관리 요구 (CEP 정지 경고)",
    },
    "[N-]=[N+]=NCc1ccc(-c2ccccc2-c2nn[nH]n2)cc1": {
        "name": "AZBT — 5-(4'-(azidomethyl)biphenyl-2-yl)-1H-tetrazole",
        "cas": "152708-24-2",
        "source": "사르탄 아지도 불순물 회수 사례에서 정량 대상으로 사용되는 표준물질",
    },
}


def main() -> int:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")

    from engines.p36_suzuki_impurity_rules import get_suzuki_predictor
    from p31_process_reaction_rules import (
        PROCESS_REACTION_RULES,
        _mols_from_smiles_list,
        _run_binary_smirks,
    )
    from p51_prime_backend_core import HybridPropertyEngine
    import p73_applicability_domain as ad_mod

    try:
        from p60_m7_empirical_db import get_default_db_connection
        conn = get_default_db_connection()
    except Exception as exc:  # noqa: BLE001
        say(f"(실측 DB 연결 실패: {exc} — Layer A 없이 진행)")
        conn = None
    from p65_m7_ensemble_matrix import run_m7_ensemble_matrix
    import p63_deepchem_predictor as dp

    say(BAR)
    say("  PRIME-DS  Case Study (DS) — 로사르탄 공정 유래 불순물 예측 및 ICH M7 판정")
    say(BAR)
    say()
    say(f"  Layer C 모델 : {dp.model_status()['model_label']}")
    say(f"  적용범위     : {ad_mod.status()['fingerprint']} · 학습셋 {ad_mod.status()['n_train']:,}건")
    say()

    rows: list[dict] = []

    def assess(step: str, rule: str, role: str, smi: str) -> dict:
        mol = Chem.MolFromSmiles(smi)
        canon = Chem.MolToSmiles(mol)
        props = HybridPropertyEngine.calculate_properties(canon) or {}
        res = run_m7_ensemble_matrix(canon, db_connection=conn)
        ad = ad_mod.applicability_domain(canon)
        rep = REPORTED.get(canon)
        alerts = "; ".join(res.layer_b.alerts) if getattr(res.layer_b, "alerts", None) else ""
        row = {
            "step": step,
            "rule": rule,
            "role": role,
            "smiles": canon,
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "mw": round(Descriptors.MolWt(mol), 2),
            "logp": props.get("LogP_Baseline", ""),
            "bp_pred_c": props.get("Boiling_Point_C", ""),
            "layer_a_hit": bool(getattr(res.layer_a, "found", False)),
            "layer_b_alerts": alerts,
            "layer_c_prob": round(float(res.layer_c.p_mutagenic), 3),
            "final_class": int(res.final_class),
            "ad_dmax": round(ad.d_max, 3) if ad.d_max is not None else "",
            "ad_status": ad.label,
            "reported_name": rep["name"] if rep else "",
            "reported_cas": rep["cas"] if rep else "",
            "reported_source": rep["source"] if rep else "",
        }
        rows.append(row)
        return row

    def show(r: dict) -> None:
        tag = "★ 문헌 보고" if r["reported_name"] else ""
        say(f"    [{r['role']:<13}] {r['rule']}")
        say(f"        {r['smiles']}")
        say(f"        {r['formula']}  MW {r['mw']:>7.2f}  LogP {r['logp']:>5}  "
            f"bp(추정) {r['bp_pred_c']:>6} °C")
        say(f"        ICH M7 Class {r['final_class']}  ·  P(Ames) {r['layer_c_prob']:.2f}  ·  "
            f"d_max {r['ad_dmax']} ({r['ad_status']})")
        if r["layer_b_alerts"]:
            say(f"        Layer B 경보 : {r['layer_b_alerts']}")
        if tag:
            say(f"        {tag} — {r['reported_name']} (CAS {r['reported_cas']})")
            say(f"                  {r['reported_source']}")
        say()

    # ── 1단계 ────────────────────────────────────────────────────────────
    say(BAR)
    say("  1단계 · Suzuki-Miyaura 커플링 — OTBN (2-cyano-4'-methylbiphenyl) 생성")
    say(BAR)
    say(f"    아릴 할라이드 : {SUZUKI_HALIDE[0]}  {SUZUKI_HALIDE[1]}")
    say(f"    보론산        : {SUZUKI_BORONIC[0]}  {SUZUKI_BORONIC[1]}")
    say()

    bundle = get_suzuki_predictor().predict_all_scenarios_marshaled(
        SUZUKI_HALIDE[1], SUZUKI_BORONIC[1]
    )
    for pred in bundle.get("predicted_impurities") or []:
        role = "Main product" if pred.get("type") == "Main" else "Side impurity"
        show(assess("1 · Suzuki coupling", pred.get("reaction_name", ""), role, pred["smiles"]))

    # ── 2단계 ────────────────────────────────────────────────────────────
    say(BAR)
    say("  2단계 · 테트라졸 형성 단계의 아자이드 친핵 치환")
    say(BAR)
    say(f"    아자이드원 : {AZIDE[0]}  {AZIDE[1]}")
    say()

    rule = PROCESS_REACTION_RULES["Nucleophilic_Substitution_with_Azide"]
    for sub_name, sub_smi in AZIDE_SUBSTRATES:
        say(f"    기질 : {sub_name}")
        say(f"           {sub_smi}")
        say()
        mols = _mols_from_smiles_list([sub_smi, AZIDE[1]])
        if len(mols) < 2:
            say("      (반응물 해석 실패)")
            continue
        ma, mb = mols[0], mols[1]
        mains = _run_binary_smirks(ma, mb, rule["main_smirks"])
        sides = _run_binary_smirks(ma, mb, rule["side_smirks"])
        for smi in mains:
            show(assess("2 · Azide substitution", f"Azide substitution on {sub_name}",
                        "Process impurity", smi))
        for smi in sides:
            show(assess("2 · Azide substitution", f"Competing E2 elimination on {sub_name}",
                        "Side impurity", smi))

    # ── 대조군 — 로사르탄 자체 ───────────────────────────────────────────
    say(BAR)
    say("  대조 · 로사르탄 원약 자체")
    say(BAR)
    say()
    show(assess("0 · Reference", "Drug substance (reference)", "API", LOSARTAN))

    # ── 요약 ─────────────────────────────────────────────────────────────
    say(BAR)
    say("  요약")
    say(BAR)
    gen = [r for r in rows if r["step"] != "0 · Reference"]
    imp = [r for r in gen if r["role"] != "Main product"]
    c3 = [r for r in imp if r["final_class"] == 3]
    matched = [r for r in gen if r["reported_name"]]
    say(f"    생성 구조 (주산물 포함) : {len(gen)}종")
    say(f"    불순물                  : {len(imp)}종")
    say(f"    ICH M7 Class 3 (변이원성 우려) : {len(c3)}종")
    for r in c3:
        say(f"        · {r['reported_name'] or r['smiles']}  ← {r['layer_b_alerts']}")
    say(f"    문헌·규제 문서 보고 불순물과 일치 : {len(matched)}종")
    for r in matched:
        say(f"        · {r['reported_name']} (CAS {r['reported_cas']})  Class {r['final_class']}")
    say()

    DATA_DIR.mkdir(exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    OUT_TXT.write_text("\n".join(_log), encoding="utf-8")
    print(f"저장: {OUT_CSV}")
    print(f"저장: {OUT_TXT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
