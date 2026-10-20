"""
논문 Phase 4 — 생성 이득 정량화.

네 개 사례연구에서 플랫폼이 «사람 대신 떠올려 그려 준» 구조가 몇 개이고,
그중 몇 개가 ICH M7 상 관리 대상이며, 몇 개가 문헌·규제 문서에 실제로
보고된 것인지 센다.

핵심 셈법
---------
Derek/Sarah 같은 구조 입력형 시스템은 «주어진 구조»를 판정할 뿐 구조를
제안하지 못한다. 따라서 같은 결과를 얻으려면 아래 «고유 구조» 수만큼을
사람이 직접 떠올려 그려 넣어야 한다. 그 수가 곧 선제적 생성의 이득이다.

결과물
------
    data/generation_benefit_report.txt
    data/generation_benefit_table.csv
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

import os

# ── 경로 ──────────────────────────────────────────────────────────────
#  BASE  : 이 스크립트가 있는 폴더 = 프로젝트 루트
#  PAPER : 논문자료 폴더. 다른 곳에 두셨다면 환경변수로 알려 주십시오.
#            set PRIME_PAPER_DIR=<논문자료 폴더의 전체 경로>
BASE = Path(__file__).resolve().parent
PAPER = Path(os.environ.get("PRIME_PAPER_DIR") or (BASE / "논문자료"))
SRC_V = PAPER / "02_검증결과"      # 검증 산출물
SRC_C = PAPER / "03_사례연구"      # 사례연구 산출물
SRC_M = PAPER / "05_원고"          # 원고·그림·보충자료


def _need(path: Path) -> Path:
    """입력이 없으면 어디를 찾았는지 밝히고 멈춘다."""
    if not path.exists():
        raise SystemExit(
            f"\n  [실패] 입력 파일이 없습니다: {path}\n"
            f"         논문자료 폴더가 다른 곳이면 PRIME_PAPER_DIR 로 지정하십시오.\n"
        )
    return path

# 산출물은 다른 사례연구 결과와 같은 곳에 쓴다
OUT_TXT = SRC_C / "generation_benefit_report.txt"
OUT_CSV = SRC_C / "generation_benefit_table.csv"

BAR = "=" * 78
_log: list[str] = []


def say(m: str = "") -> None:
    print(m)
    _log.append(m)


def keys(series) -> set[str]:
    out = set()
    for x in series:
        if not isinstance(x, str) or not x.strip():
            continue
        m = Chem.MolFromSmiles(x)
        if m is not None:
            out.add(Chem.MolToInchiKey(m))
    return out


# 상호작용 규칙 총수 — 원고 §2.5 서술과 같은 값
N_INTERACTION_RULES = 14

# ── 자료 읽기 ───────────────────────────────────────────────────────────
ds = pd.read_csv(_need(SRC_C / "ds_process_case_study_results.csv"))
cs = pd.read_csv(_need(SRC_C / "case_study_results.csv"))
ex = pd.read_csv(_need(SRC_C / "excipient_case_study_results.csv"))
bi = pd.read_csv(_need(SRC_C / "binary_compatibility_pairs.csv"))

ds_imp = ds[(ds["step"] != "0 · Reference") & (ds["role"] != "Main product")]
cs_p = cs[cs["product_smiles"].notna()]
bi_p = bi[bi["Impurity_SMILES"].notna()]

k1 = keys(ds_imp["smiles"])
k2 = keys(cs_p["product_smiles"])
k3 = keys(ex["smiles"])
k4 = keys(bi_p["Impurity_SMILES"])

# ── 관리 대상 (Class 1 또는 3) ──────────────────────────────────────────
c1_conc = keys(ds_imp[ds_imp["final_class"].astype(int) <= 3]["smiles"])
c2_conc = keys(cs_p[cs_p["final_class"].isin(["Class 1", "Class 3"])]["product_smiles"])
c3_conc = keys(ex[ex["ich_m7"].isin(["Class 1", "Class 3"])]["smiles"])

# ── 문헌·규제 문서에 보고된 것 ──────────────────────────────────────────
DOCUMENTED = {
    "2,2'-dicyanobiphenyl": ("N#Cc1ccccc1-c1ccccc1C#N", "US 6,153,810"),
    "4,4'-dimethylbiphenyl": ("Cc1ccc(-c2ccc(C)cc2)cc1", "US 6,392,080"),
    "losartan azido impurity": (
        "CCCCc1nc(Cl)c(CN=[N+]=[N-])n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
        "EDQM 2021 (Ames 양성 · TTC 관리)"),
    "azidomethyl biphenyl tetrazole": (
        "[N-]=[N+]=NCc1ccc(-c2ccccc2-c2nn[nH]n2)cc1",
        "사르탄 아지도 불순물 표준물질"),
    "NDMA": ("CN(C)N=O", "FDA 회수 — 라니티딘 · 니자티딘 (2019–2020)"),
}
doc_keys = {}
for name, (smi, src) in DOCUMENTED.items():
    doc_keys[Chem.MolToInchiKey(Chem.MolFromSmiles(smi))] = (name, src)

allk = k1 | k2 | k3 | k4
conc = c1_conc | c2_conc | c3_conc
doc_hit = {k: v for k, v in doc_keys.items() if k in allk}

# ── 표 ──────────────────────────────────────────────────────────────────
ROWS = [
    {
        "case": "1 — DS process (losartan route)",
        "user_supplies": "3 starting materials, 2 steps",
        "rule_applications": 5 + 2 * 2,
        "records": len(ds_imp),
        "distinct": len(k1),
        "controlled": len(c1_conc),
        "documented": len([k for k in k1 if k in doc_keys]),
    },
    {
        "case": "2 — NDSRI formation",
        "user_supplies": "4 APIs + nitrite source",
        "rule_applications": cs["api_name"].nunique(),
        "records": len(cs_p),
        "distinct": len(k2),
        "controlled": len(c2_conc),
        "documented": len([k for k in k2 if k in doc_keys]),
    },
    {
        "case": "3 — excipient matrix",
        "user_supplies": "3 APIs × 5 excipients",
        "rule_applications": 3 * 5 * N_INTERACTION_RULES,
        "records": len(ex),
        "distinct": len(k3),
        "controlled": len(c3_conc),
        "documented": len([k for k in k3 if k in doc_keys]),
    },
    {
        "case": "4 — binary 1:1 pairs",
        "user_supplies": "15 API–excipient pairs",
        "rule_applications": 15 * N_INTERACTION_RULES,
        "records": len(bi_p),
        "distinct": len(k4),
        "controlled": None,          # 이원 연구는 위험 점수 체계라 M7 등급을 매기지 않는다
        "documented": len([k for k in k4 if k in doc_keys]),
    },
]

say(BAR)
say("  PRIME-DS  Phase 4 — 생성 이득 정량화")
say(BAR)
say()
say("  구조 입력형 시스템(Derek/Sarah)은 주어진 구조를 판정할 뿐 구조를 제안하지")
say("  못한다. 아래 «고유 구조» 열이 곧 같은 결과를 얻기 위해 사람이 직접")
say("  떠올려 그려 넣어야 하는 구조의 수다.")
say()
say(f"  {'사례':<34}{'규칙 적용':>10}{'산출 레코드':>12}{'고유 구조':>10}{'관리 대상':>10}{'문헌 보고':>10}")
say("  " + "-" * 88)
for r in ROWS:
    ctrl = "—" if r["controlled"] is None else f"{r['controlled']}"
    say(f"  {r['case']:<34}{r['rule_applications']:>10}{r['records']:>12}"
        f"{r['distinct']:>10}{ctrl:>10}{r['documented']:>10}")
say("  " + "-" * 88)
say(f"  {'합계 (중복 제거)':<34}{sum(r['rule_applications'] for r in ROWS):>10}"
    f"{sum(r['records'] for r in ROWS):>12}{len(allk):>10}{len(conc):>10}{len(doc_hit):>10}")
say()

say(f"  · 규칙 적용 {sum(r['rule_applications'] for r in ROWS)}회 → 고유 구조 {len(allk)}개")
say(f"  · 고유 구조 {len(allk)}개 중 ICH M7 Class 1 또는 3 : {len(conc)}개 "
    f"({100 * len(conc) / len(allk):.0f}%)")
say(f"  · 관리 대상 1건을 찾기 위해 검토한 구조 : {len(allk) / len(conc):.1f}개")
say(f"  · 문헌·규제 문서에 실제로 보고된 구조 : {len(doc_hit)}개")
for k, (name, src) in doc_hit.items():
    say(f"        · {name} — {src}")
say()

say("  사례 3 과 4 의 중복 — 같은 화학을 두 방식으로 평가한 것이므로 구조가 겹친다.")
say(f"        사례 3 고유 {len(k3)} · 사례 4 고유 {len(k4)} · 합집합 {len(k3 | k4)} "
    f"(겹침 {len(k3 & k4)})")
say(f"  사례 2 의 4 레코드는 고유 구조 {len(k2)} 개다 — NDMA 가 라니티딘·니자티딘")
say("        양쪽에서 나오기 때문이며, 이것이 두 품목이 같은 이유로 회수된 까닭이다.")
say()
say(BAR)

SRC_C.mkdir(parents=True, exist_ok=True)
with OUT_CSV.open("w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(ROWS[0].keys()))
    w.writeheader()
    w.writerows(ROWS)
OUT_TXT.write_text("\n".join(_log), encoding="utf-8")
print(f"\n저장: {OUT_CSV}")
print(f"저장: {OUT_TXT}")
