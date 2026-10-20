"""
오픈소스 부형제 레지스트리 — FDA IID(비활성 성분) + PubChem 구조 연동.

- 로컬 마스터(USP/NF, Handbook) + FDA IID 보충 목록 병합
- SMILES 미보유 항목은 PubChem REST API로 lazy 해석 (캐시)
"""

# 목록을 76종으로 늘림.
# FDA Inactive Ingredient Database에서 온 29종을 더함. 구조가 없는 항목은 PubChem에서 받아옴(resolve_pubchem_smiles).

from __future__ import annotations

import logging
from functools import lru_cache
from typing import Any

from p70_excipient_master_library import COMPLETE_EXCIPIENT_MASTER, EXCIPIENT_CATEGORY_ORDER

_LOGGER = logging.getLogger(__name__)

# FDA Inactive Ingredient Database (IID) · USP-NF Supplement — PubChem명 매핑
# 출처: https://www.accessdata.fda.gov/scripts/cder/iig/index.htm
#
# ===========================================================================
# ⚠️ 2026-09-03 — 30종 중 25종이 `smiles`·`reactive_probes` 없이 등록되어
#    있었습니다. 그 결과 `probe_smiles_for_entry`가 **매트릭스 평가 루프
#    안에서** `resolve_pubchem_smiles`(네트워크)를 호출했습니다.
#
#      · 첫 평가 한 번에 최대 25회 외부 API 왕복
#      · 실패하면 `except Exception: return None` 으로 **무음** 처리되어
#        해당 부형제가 조용히 평가에서 빠짐 (위음성)
#      · 오프라인에서는 재현 자체가 불가능
#
#    이제 **전 항목의 구조를 소스에 직접 채워** 넣었습니다. 반응 예측
#    경로에서 네트워크를 타지 않습니다. PubChem 조회는 소스에 구조가 없는
#    항목에 대한 마지막 수단으로만 남고, 실패는 무음이 아니라 로그로 남습니다.
#
#    구조는 p19의 STRUCTURE_PATCH와 같은 값이며 전부 RDKit 파싱을
#    검증했습니다. `mw_kind`는 p70과 같은 세 값 중 하나입니다
#    ("molecular" · "monomer_surrogate" · "undefined").
#    고분자·천연 혼합물의 `smiles`는 **대표 반복 단위/조각**이며
#    몰비 계산에 쓰지 않습니다.
# ===========================================================================
_OPENSOURCE_SUPPLEMENT: dict[str, dict[str, Any]] = {
    "Gelatin": {
        "pubchem_name": "gelatin",
        # 대표 반복 단위 (Gly-Pro-Hyp) — 유리 아미노기가 반응 자리
        "smiles": "NCC(=O)N1CCCC1C(=O)NCC(=O)O",
        "reactive_probes": ["NCC(=O)O", "NCCCCC(N)C(=O)O"],
        "mw_kind": "monomer_surrogate",
        "structure_note": "가수분해 콜라겐 단백질 — 등재 구조는 Gly-Pro-Hyp 삼펩타이드 대용입니다.",
        "category": "Binder/Diluent",
        "label_ko": "젤라틴 (Gelatin)",
        "alerts": "단백질 가수분해 — 아민기 API Maillard · 이성질화 리스크",
        "source": "FDA IID / PubChem",
    },
    "Xanthan_Gum": {
        "pubchem_name": "xanthan gum",
        "smiles": "OCC1OC(O)C(O)C(O)C1O",
        "reactive_probes": ["OC(=O)C1OC(O)C(O)C(O)C1O"],
        "mw_kind": "monomer_surrogate",
        "structure_note": "다당류 고분자 — 등재 구조는 글루코스 반복 단위 대용입니다.",
        "category": "Binder/Diluent",
        "label_ko": "잔탄검 (Xanthan gum)",
        "alerts": "다당류 — 수분 보유 · 고체상 가수분해 촉진",
        "source": "FDA IID / PubChem",
    },
    "Carboxymethylcellulose_Sodium": {
        "pubchem_name": "carboxymethyl cellulose sodium",
        # 대이온 [Na+] 포함 — 순전하 0 (C8H13NaO8)
        "smiles": "[Na+].OCC1OC(O)C(OCC(=O)[O-])C(O)C1O",
        "reactive_probes": ["OCC(=O)O"],
        "mw_kind": "monomer_surrogate",
        "structure_note": "카복시메틸셀룰로오스 나트륨 — 등재 구조는 CM-글루코스 반복 단위입니다.",
        "category": "Binder/Diluent",
        "label_ko": "카르복시메틸셀룰로오스 나트륨 (CMC)",
        "alerts": "점증제 — 국소 수분 활성도 증가",
        "source": "FDA IID / PubChem",
    },
    "Alginic_Acid": {
        "pubchem_name": "alginic acid",
        "smiles": "OC1C(O)C(OC(C1O)C(=O)O)O",
        "reactive_probes": ["OC1C(O)C(O)OC(C(=O)O)C1O"],
        "mw_kind": "monomer_surrogate",
        "structure_note": "만누론산/굴루론산 공중합 다당류 — 등재 구조는 우론산 반복 단위입니다.",
        "category": "Binder/Diluent",
        "label_ko": "알긴산 (Alginic acid)",
        "alerts": "다당류 산 — pH 완충 · 이온 교환",
        "source": "FDA IID / PubChem",
    },
    "Acacia_Gum": {
        "pubchem_name": "gum arabic",
        "smiles": "OCC1OC(O)C(O)C(O)C1O",
        "reactive_probes": ["O=CC(O)C(O)C(O)CO"],
        "mw_kind": "monomer_surrogate",
        "structure_note": "아라비노갈락탄 단백질 복합체 — 등재 구조는 당 반복 단위 대용입니다.",
        "category": "Binder/Diluent",
        "label_ko": "아라비아검 (Acacia / Gum arabic)",
        "alerts": "천연 다당류 — Maillard 보조",
        "source": "FDA IID / PubChem",
    },
    "Sodium_Chloride": {
        "pubchem_name": "sodium chloride",
        "smiles": "[Na+].[Cl-]",
        "reactive_probes": [],
        "mw_kind": "molecular",
        "category": "Binder/Diluent",
        "label_ko": "염화나트륨 (NaCl)",
        "alerts": "이온강도 변화 — 용해도·가수분해 속도론 영향",
        "source": "FDA IID / PubChem",
    },
    "Dextrose": {
        "pubchem_name": "D-glucose",
        "smiles": "OCC1OC(O)C(O)C(O)C1O",
        "reactive_probes": ["O=CC(O)CO"],
        "mw_kind": "molecular",
        "category": "Binder/Diluent",
        "label_ko": "덱스트로스 (D-Glucose)",
        "alerts": "🔥 환원당 — Maillard 불순물",
        "source": "FDA IID / PubChem",
    },
    "Calcium_Phosphate_Tribasic": {
        "pubchem_name": "tricalcium phosphate",
        "smiles": "[Ca+2].[Ca+2].[Ca+2].[O-]P([O-])([O-])=O.[O-]P([O-])([O-])=O",
        "reactive_probes": [],
        "mw_kind": "molecular",
        "category": "Binder/Diluent",
        "label_ko": "인산칼슘 (Tribasic calcium phosphate)",
        "alerts": "알칼리성 희석제 — pH 상승",
        "source": "FDA IID / PubChem",
    },
    "Titanium_Dioxide": {
        "pubchem_name": "titanium dioxide",
        "smiles": "O=[Ti]=O",
        "reactive_probes": [],
        "mw_kind": "undefined",
        "structure_note": "TiO₂ 망상 고체 — 분자량 정의 없음.",
        "category": "Binder/Diluent",
        "label_ko": "이산화티타늄 (TiO₂)",
        "alerts": "광촉매 표면 — 광분해 보조 (제형 코팅)",
        "source": "FDA IID / PubChem",
    },
    "Glyceryl_Monostearate": {
        "pubchem_name": "glyceryl monostearate",
        "smiles": "CCCCCCCCCCCCCCCCCC(=O)OCC(O)CO",
        # 에스터 → 에스터 교환 · 아미놀리시스, 가수분해 시 스테아르산
        "reactive_probes": ["CCCCCCCCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCCCCCC(=O)O"],
        "mw_kind": "molecular",
        "category": "Lubricant/Glidant",
        "label_ko": "글리세릴 모노스테아레이트",
        "alerts": "지방산 에스터 — 가수분해 · 산화",
        "source": "FDA IID / PubChem",
    },
    "Hydrogenated_Castor_Oil": {
        "pubchem_name": "hydrogenated castor oil",
        "smiles": "CCCCCCC(O)CCCCCCCCCCC(=O)OCC(O)CO",
        "reactive_probes": ["CCCCCCC(O)CCCCCCCCCCC(=O)O"],
        "mw_kind": "undefined",
        "structure_note": "수소첨가 트리글리세리드 혼합물 — 등재 구조는 12-히드록시스테아르산 모노에스터 대용입니다.",
        "category": "Lubricant/Glidant",
        "label_ko": "수소첨가 피마자유 (Castor oil, hydrogenated)",
        "alerts": "유지 — 산화·과산화물",
        "source": "FDA IID / PubChem",
    },
    "Light_Mineral_Oil": {
        "pubchem_name": "mineral oil",
        "smiles": "CCCCCCCCCCCCCCCC",
        "reactive_probes": [],
        "mw_kind": "undefined",
        "structure_note": "포화 탄화수소 혼합물 — 등재 구조는 헥사데칸 대표 사슬입니다.",
        "category": "Lubricant/Glidant",
        "label_ko": "미네랄오일 (Light mineral oil)",
        "alerts": "탄화수소 — 추출·용매 잔류",
        "source": "FDA IID / PubChem",
    },
    "Sodium_Benzoate": {
        "pubchem_name": "sodium benzoate",
        "smiles": "[Na+].[O-]C(=O)c1ccccc1",
        "reactive_probes": ["OC(=O)c1ccccc1"],
        "mw_kind": "molecular",
        "category": "Preservative/Antioxidant",
        "label_ko": "벤조산나트륨",
        "alerts": "방부제 — 산성 pH · 벤조산 유리",
        "source": "FDA IID / PubChem",
    },
    "Potassium_Sorbate": {
        "pubchem_name": "potassium sorbate",
        # α,β-불포화 → Michael 수용체
        "smiles": "[K+].C/C=C/C=C/C(=O)[O-]",
        "reactive_probes": ["C/C=C/C=C/C(=O)O"],
        "mw_kind": "molecular",
        "category": "Preservative/Antioxidant",
        "label_ko": "소르빈산칼륨",
        "alerts": "불포화 산 염 — 산화 분해",
        "source": "FDA IID / PubChem",
    },
    "Butylated_Hydroxyanisole_BHA": {
        "pubchem_name": "butylated hydroxyanisole",
        "smiles": "COc1ccc(O)c(C(C)(C)C)c1",
        "reactive_probes": ["COc1ccc(O)c(C(C)(C)C)c1"],
        "mw_kind": "molecular",
        "category": "Preservative/Antioxidant",
        "label_ko": "BHA (Butylated hydroxyanisole)",
        "alerts": "항산화제 — quinone 부산물",
        "source": "FDA IID / PubChem",
    },
    "Sodium_Metabisulfite": {
        "pubchem_name": "sodium metabisulfite",
        # ⚠️ Na2S2O5 — S–O–S. S–S 로 쓰면 dithionite(Na2S2O4)가 됩니다.
        "smiles": "[Na+].[Na+].[O-]S(=O)OS(=O)[O-]",
        "reactive_probes": ["OS(=O)O", "[O-]S(=O)O"],
        "mw_kind": "molecular",
        "category": "Preservative/Antioxidant",
        "label_ko": "피메타중아황산나트륨",
        "alerts": "🔥 환원제 — 아민·알데하이드 API 반응",
        "source": "FDA IID / PubChem",
    },
    "Hypromellose_Phthalate": {
        "pubchem_name": "hypromellose phthalate",
        "smiles": "OC(=O)c1ccccc1C(=O)OCC1OC(O)C(O)C(O)C1O",
        "reactive_probes": ["OC(=O)c1ccccc1C(=O)O", "OC(=O)c1ccccc1C(=O)OC"],
        "mw_kind": "monomer_surrogate",
        "structure_note": "셀룰로오스 에터 프탈산 에스터 — 등재 구조는 프탈로일 글루코스 반복 단위입니다.",
        "category": "Polymer/Film Coating",
        "label_ko": "프탈산히프로멜로오스 (HPMCP)",
        "alerts": "장용 코팅 — pH 의존 용출",
        "source": "FDA IID / PubChem",
    },
    "Cellulose_Acetate_Phthalate": {
        "pubchem_name": "cellulose acetate phthalate",
        # 가수분해로 아세트산 · 프탈산 유리
        "smiles": "CC(=O)OCC1OC(O)C(OC(=O)c2ccccc2C(=O)O)C(O)C1O",
        "reactive_probes": ["CC(=O)O", "OC(=O)c1ccccc1C(=O)O", "CC(=O)OC"],
        "mw_kind": "monomer_surrogate",
        "structure_note": "셀룰로오스 아세테이트 프탈레이트 — 등재 구조는 대표 반복 단위입니다.",
        "category": "Polymer/Film Coating",
        "label_ko": "프탈산셀룰로오스아세테이트 (CAP)",
        "alerts": "장용 코팅 — 가수분해 잔류 아세트산",
        "source": "FDA IID / PubChem",
    },
    "Povidone_K90": {
        "pubchem_name": "povidone K90",
        # 대표 반복 단위(N-비닐피롤리돈). K30과 같은 값입니다.
        "smiles": "CCN1CCCC1=O",
        "reactive_probes": ["C(=O)O", "OO"],
        "mw_kind": "monomer_surrogate",
        "structure_note": "선형 PVP 고분자. K값은 분자량 등급이며 반복 단위는 등급과 무관하게 같습니다.",
        "category": "Polymer/Film Coating",
        "label_ko": "포비돈 K90",
        "alerts": "고분자 결합제 — 과산화물 잔류",
        "source": "FDA IID / PubChem",
    },
    "Polysorbate_60": {
        "pubchem_name": "polysorbate 60",
        "smiles": "CCCCCCCCCCCCCCCCCC(=O)OCCOCCO",
        "reactive_probes": ["C(=O)O", "OO"],
        "mw_kind": "monomer_surrogate",
        "structure_note": "폴리옥시에틸렌 소르비탄 모노스테아레이트 — 등재 구조는 대표 조각입니다.",
        "category": "Surfactant",
        "label_ko": "폴리소르베이트 60 (Tween 60)",
        "alerts": "에스터 계면활성제 — 과산화물",
        "source": "FDA IID / PubChem",
    },
    "Lecithin": {
        "pubchem_name": "lecithin",
        # 포스파티딜콜린 — 에스터 결합 + 4급 암모늄 (양쪽성 이온, 순전하 0)
        "smiles": "CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC",
        "reactive_probes": ["CCCCCCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCCCC(=O)O"],
        "mw_kind": "undefined",
        "structure_note": "인지질 혼합물 — 등재 구조는 디팔미토일 포스파티딜콜린 대표입니다.",
        "category": "Surfactant",
        "label_ko": "레시틴 (Lecithin)",
        "alerts": "인지질 — 산화·라디칼",
        "source": "FDA IID / PubChem",
    },
    "Sodium_Phosphate_Dibasic": {
        "pubchem_name": "disodium hydrogen phosphate",
        "smiles": "[Na+].[Na+].OP([O-])([O-])=O",
        "reactive_probes": [],
        "mw_kind": "molecular",
        "category": "Buffer/pH Modifier",
        "label_ko": "인산수소나트륨",
        "alerts": "완충염 — pH 7~8 염기성 가수분해",
        "source": "FDA IID / PubChem",
    },
    "Potassium_Phosphate": {
        "pubchem_name": "potassium phosphate",
        "smiles": "[K+].[K+].OP([O-])([O-])=O",
        "reactive_probes": [],
        "mw_kind": "molecular",
        "category": "Buffer/pH Modifier",
        "label_ko": "인산칼륨",
        "alerts": "완충 — 이온강도·pH 변동",
        "source": "FDA IID / PubChem",
    },
    "Sodium_Hydroxide": {
        "pubchem_name": "sodium hydroxide",
        "smiles": "[Na+].[OH-]",
        "reactive_probes": [],
        "mw_kind": "molecular",
        "category": "Buffer/pH Modifier",
        "label_ko": "수산화나트륨 (NaOH)",
        "alerts": "강염기 — 염기성 가수분해 촉진",
        "source": "FDA IID / PubChem",
    },
    "Hydrochloric_Acid": {
        "pubchem_name": "hydrochloric acid",
        "smiles": "Cl",
        "reactive_probes": [],
        "mw_kind": "molecular",
        "category": "Buffer/pH Modifier",
        "label_ko": "염산 (HCl)",
        "alerts": "강산 — 산 가수분해 촉진",
        "source": "FDA IID / PubChem",
    },
    "Saccharin_Sodium": {
        "pubchem_name": "saccharin sodium",
        "smiles": "[Na+].O=C1[N-]S(=O)(=O)c2ccccc21",
        "reactive_probes": [],
        "mw_kind": "molecular",
        "category": "Sweetener/Flavor",
        "label_ko": "사카린나트륨",
        "alerts": "유기산 염 — pH·이온쌍 형성",
        "source": "FDA IID / PubChem",
    },
    "Peppermint_Oil": {
        "pubchem_name": "menthol",
        # 대표 성분 멘톨 + 멘톤(케톤 → 중아황산 부가 가능)
        "smiles": "CC(C)C1CCC(C)CC1O",
        "reactive_probes": ["CC(C)C1CCC(C)CC1=O"],
        "mw_kind": "undefined",
        "structure_note": "테르펜 혼합 정유 — 등재 구조는 대표 성분 멘톨입니다.",
        "category": "Sweetener/Flavor",
        "label_ko": "페퍼민트오일 (대표: 멘톨)",
        "alerts": "terpene — 산화 · API 추출",
        "source": "FDA IID / PubChem",
    },
    "Isopropyl_Alcohol": {
        "pubchem_name": "isopropyl alcohol",
        "smiles": "CC(C)O",
        # ⚠️ IPA는 2급 알코올 → 산화 산물은 **아세톤**입니다.
        #    아세트알데하이드(CC=O)는 에탄올 유래이므로 여기 쓰면 위양성입니다.
        "reactive_probes": ["CC(C)=O"],
        "mw_kind": "molecular",
        "category": "Solvent/Co-solvent",
        "label_ko": "이소프로필알코올 (IPA)",
        "alerts": "용매 — 산화 → 아세톤",
        "source": "FDA IID / PubChem",
    },
    "Benzyl_Benzoate": {
        "pubchem_name": "benzyl benzoate",
        "smiles": "O=C(OCc1ccccc1)c1ccccc1",
        "reactive_probes": ["O=C(OC)c1ccccc1", "O=C(O)c1ccccc1"],
        "mw_kind": "molecular",
        "category": "Solvent/Co-solvent",
        "label_ko": "벤질벤조산에스테르",
        "alerts": "에스터 용매 — 가수분해 부산물",
        "source": "FDA IID / PubChem",
    },
}


@lru_cache(maxsize=512)
def resolve_pubchem_smiles(pubchem_name: str) -> str | None:
    """PubChem 이름 → SMILES. 소스에 구조가 없을 때의 **마지막 수단**입니다.

    2026-09-03 수정: 전에는 실패를 `except Exception: return None` 으로
    조용히 삼켰습니다. 네트워크가 없으면 해당 부형제가 «이상 없음»으로
    보고되는데 사용자는 조회가 실패했다는 사실 자체를 알 수 없었습니다.
    이제 실패를 로그로 남깁니다(동작은 그대로 None 반환).
    """
    if not pubchem_name:
        return None
    try:
        from p25_pubchem_client import fetch_smiles_by_compound_name

        resolved = fetch_smiles_by_compound_name(pubchem_name)
    except Exception as exc:
        _LOGGER.warning(
            "PubChem 조회 실패 (%s): %s: %s — 이 부형제는 구조 없이 평가에서 빠집니다",
            pubchem_name, type(exc).__name__, exc,
        )
        return None
    if not resolved:
        _LOGGER.warning(
            "PubChem 조회 결과 없음 (%s) — 이 부형제는 구조 없이 평가에서 빠집니다",
            pubchem_name,
        )
        return None
    return resolved


def validate_registry_smiles(*, log: bool = True) -> dict[str, Any]:
    """병합 라이브러리 전 항목의 구조를 검증합니다 (로드 시점 점검용).

    2026-09-03 추가. 마그네슘 스테아레이트의 SMILES가 문법 오류였는데도
    `probe_smiles_for_entry`가 None을 조용히 버려 **전 계층이 무음**
    이었습니다. 이제 실패를 명시적으로 셉니다.
    """
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.error")

    lib = get_full_excipient_library()
    smiles_ok = 0
    smiles_failed: list[tuple[str, str]] = []
    probe_failed: list[tuple[str, str]] = []
    no_structure: list[str] = []

    for key, entry in lib.items():
        smi = str(entry.get("smiles") or "").strip()
        if not smi:
            no_structure.append(key)
            if log:
                _LOGGER.warning(
                    "부형제 %s: 소스에 구조가 없습니다 — 평가 중 PubChem 네트워크 조회가 발생합니다",
                    key,
                )
        elif Chem.MolFromSmiles(smi) is None:
            smiles_failed.append((key, smi))
            if log:
                _LOGGER.warning("부형제 %s: smiles 파싱 실패 %r", key, smi)
        else:
            smiles_ok += 1

        for probe in entry.get("reactive_probes") or []:
            if Chem.MolFromSmiles(str(probe)) is None:
                probe_failed.append((key, str(probe)))
                if log:
                    _LOGGER.warning("부형제 %s: reactive_probe 파싱 실패 %r", key, probe)

    return {
        "total": len(lib),
        "smiles_ok": smiles_ok,
        "smiles_failed": smiles_failed,
        "probe_failed": probe_failed,
        "no_structure": no_structure,
        "network_lookups_required": len(no_structure),
    }


def _normalize_entry(data: dict[str, Any]) -> dict[str, Any]:
    out = dict(data)
    out.setdefault("smiles", "")
    out.setdefault("reactive_probes", [])
    out["risk_alerts"] = out.get("alerts") or out.get("risk_alerts") or ""
    out.setdefault("source", "PRIME-DS Master")
    return out


def get_full_excipient_library() -> dict[str, dict[str, Any]]:
    """로컬 마스터 + FDA IID 오픈소스 보충 병합."""
    merged: dict[str, dict[str, Any]] = {}
    for key, data in COMPLETE_EXCIPIENT_MASTER.items():
        merged[key] = _normalize_entry(data)
    for key, data in _OPENSOURCE_SUPPLEMENT.items():
        if key not in merged:
            entry = _normalize_entry(data)
            entry["source"] = data.get("source", "FDA IID / PubChem")
            merged[key] = entry
    return merged


def probe_smiles_for_entry(entry: dict[str, Any]) -> list[str]:
    """부형제 엔트리 → RDKit 반응 투입용 SMILES 목록.

    2026-09-03 수정 — 두 가지가 무음이었습니다.
      ① 파싱 실패한 SMILES를 아무 말 없이 버렸습니다. 마그네슘
         스테아레이트의 문법 오류가 여기서 조용히 사라졌습니다.
      ② 소스에 구조가 없으면 **평가 루프 안에서** PubChem을 호출했습니다.
    이제 ①은 로그로 남기고, ②는 소스에 구조를 채워 사실상 발생하지
    않지만 남은 경우를 대비해 경고를 남깁니다.
    """
    from rdkit import Chem

    label = str(entry.get("label_ko") or entry.get("pubchem_name") or "?")
    smis: list[str] = []
    seen: set[str] = set()

    def _add(smi: str, *, kind: str) -> None:
        s = (smi or "").strip()
        if not s or s in seen:
            return
        if Chem.MolFromSmiles(s) is None:
            _LOGGER.warning(
                "부형제 %s: %s 파싱 실패로 반응 평가에서 제외됩니다 — %r", label, kind, s
            )
            return
        seen.add(s)
        smis.append(s)

    base = str(entry.get("smiles") or "").strip()
    if base:
        _add(base, kind="smiles")
    for probe in entry.get("reactive_probes") or []:
        _add(str(probe), kind="reactive_probe")
    if not smis:
        pc = str(entry.get("pubchem_name") or entry.get("label_ko") or "").strip()
        if pc:
            _LOGGER.warning(
                "부형제 %s: 소스에 유효한 구조가 없어 PubChem 네트워크 조회로 넘어갑니다 (%s)",
                label, pc,
            )
            resolved = resolve_pubchem_smiles(pc)
            if resolved:
                _add(resolved, kind="pubchem 해석 결과")
    return smis


def opensource_catalog_stats() -> dict[str, int]:
    lib = get_full_excipient_library()
    master_n = len(COMPLETE_EXCIPIENT_MASTER)
    return {
        "total": len(lib),
        "master": master_n,
        "opensource": len(lib) - master_n,
        "categories": len(EXCIPIENT_CATEGORY_ORDER),
    }
