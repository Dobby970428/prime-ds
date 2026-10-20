"""
PRIME-DS Backend Core Suite — hybrid physicochemical, dual-layer process parsing,
thermodynamic degradation tree, and Hansen/ISSCAN read-across.
"""

# | `HybridPropertyEngine` | pKa · LogD · 끓는점 ⚠️ | `p50` |
# | `DualLayerReactionParser` | 공정 반응 주/부 분리 | `p91` |
# | `DegradationTreeEngine` | 분해 트리 (JSON용) | `p92` |
# | `ReadAcrossEngine` | 유사 물질 찾기 | `p96` Flask 앱 |

# HybridPropertyEngine: 분자 구조를 보고 계산하는 게 아니라 세 갈래로 나눠 정해진 숫자를 대입함.
# 분류일 뿐 예측이 아님.

# p50 → 결과 딕셔너리
#    ├─ 화면 "물리화학적 물성" 표
#    ├─ p80 규제 보고서
#    └─ p92 Backend Core HTML
# ```

from __future__ import annotations

import copy
import math
from typing import Any

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

from p23_chem_search_engine import safe_canonical_smiles
from p44_degradation_shelf_life_filter import (
    ShelfLifeChildFilter,
    is_spin_stabilized_neutral_mol,
    is_spin_stabilized_neutral_smiles,
    should_cull_metabolic_aromatic_hydroxylation,
)
from p84_structure_images import smiles_to_base64_img

# Fallback when SQLite empirical index is unavailable
HANSEN_MUTAGENICITY_DB: list[dict[str, str]] = [
    {"name": "Aniline", "smiles": "Nc1ccccc1", "ames_result": "positive", "source": "Hansen"},
    {"name": "Benzene", "smiles": "c1ccccc1", "ames_result": "negative", "source": "Hansen"},
    {"name": "p-Chloroaniline", "smiles": "Nc1ccc(Cl)cc1", "ames_result": "positive", "source": "Hansen"},
    {"name": "Phenol", "smiles": "Oc1ccccc1", "ames_result": "negative", "source": "Hansen"},
    {"name": "Nitrobenzene", "smiles": "[O-][N+](=O)c1ccccc1", "ames_result": "positive", "source": "Hansen"},
    {"name": "Formaldehyde", "smiles": "C=O", "ames_result": "positive", "source": "Hansen"},
    {"name": "Dimethyl sulfoxide", "smiles": "CS(=O)C", "ames_result": "negative", "source": "Hansen"},
]

_ACID_SMARTS = Chem.MolFromSmarts("C(=O)[OH]")
# 염기성 질소 판정.
#   !$(NC=O)   아마이드 질소 제외 (원래 !$(C=O) 로 적혀 있었으나
#              질소 자신이 C=O 인지 보는 조건이라 아무것도 걸러내지 못했음)
#   !$(N=O)    니트로소 제외
#   !$([N+][O-]) 니트로기 제외
#   !$(Nc)     아닐린류 방향족 아민 제외 — 실제 pKa가 4~5 로 약염기
#
# 2026-09-03 수정 — 주석은 `!$(Nc)` 로 아닐린류를 제외한다고 **적어 놓고
#   SMARTS에는 그 조건이 없었습니다.** 그래서 아닐린이 지방족 아민과 같은
#   pKa 9.5 로 처리됐습니다.
#
#       아닐린 짝산의 실제 pKa ≈ 4.6  →  pH 7.4에서 거의 중성
#       따라서 LogD ≈ LogP ≈ 1.27 이어야 하는데 −0.83 이 나왔습니다.
#       **분배계수로 약 130배 오차**가 p50 물성표를 거쳐 규제 보고서까지
#       전달됐습니다. 아닐린류는 유전독성 불순물에서 가장 흔한 골격입니다.
#
#   지금은 두 갈래로 나눕니다.
#       _BASE_SMARTS          지방족(강)염기성 질소 — pKa 9.5, 아릴아민 제외
#       _AROMATIC_AMINE_SMARTS 아릴아민 — pKa 4.6 (약염기)
_BASE_SMARTS = Chem.MolFromSmarts(
    "[NX3;H2,H1,H0;!$(NC=O);!$(N=O);!$([N+][O-]);!$(N=*);!$(Nc);!$(NS(=O)=O);!$(n)]"
)
# 방향족 탄소에 직접 붙은 3가 질소 (아닐린 · p-톨루이딘 · N-메틸아닐린 등).
# 아마이드·설폰아마이드·니트로소·니트로는 제외합니다.
_AROMATIC_AMINE_SMARTS = Chem.MolFromSmarts(
    "[NX3;H2,H1,H0;$(Nc);!$(NC=O);!$(N=O);!$([N+][O-]);!$(N=*);!$(NS(=O)=O)]"
)

# pKa 대입값 (짝산 기준)
_PKA_CARBOXYLIC_ACID = 4.5
_PKA_ALIPHATIC_AMINE = 9.5
# 아닐리늄 이온의 실측 pKa. 전자끄는기가 붙으면 더 낮아지지만
# 이 엔진은 치환기 효과를 모델링하지 않으므로 모체 값을 씁니다.
_PKA_AROMATIC_AMINE = 4.6
_NITROSO_SMARTS = Chem.MolFromSmarts("[#7][#7]=[#8]")
_NITROSO_FRAGMENT_SMARTS = Chem.MolFromSmarts("[#7]=[#8]")
_MIN_DEGRADANT_MW = 40.0
READ_ACROSS_MIN_SIMILARITY = 0.5


class HybridPropertyEngine:
    """Henderson-Hasselbalch LogD + weighted BP regression (thermodynamic layer)."""

    @staticmethod
    def calculate_properties(smiles: str, *, target_ph: float = 7.4) -> dict[str, float] | None:
        mol = Chem.MolFromSmiles((smiles or "").strip())
        if mol is None:
            return None

        mw = Descriptors.MolWt(mol)
        logp_baseline = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)

        boiling_point = (195.0 + 1.25 * mw + 1.8 * tpsa + 25.0 * hbd) - 273.15

        pka = 7.0
        is_acidic = False
        is_basic = False
        is_nitrosamine = bool(_NITROSO_SMARTS and mol.HasSubstructMatch(_NITROSO_SMARTS))

        # 2026-09-03: 지방족 아민과 방향족 아민(아닐린류)을 갈라
        #   각각 pKa 9.5 / 4.6 을 적용합니다. 전에는 둘 다 9.5였습니다.
        if _ACID_SMARTS and mol.HasSubstructMatch(_ACID_SMARTS):
            pka = _PKA_CARBOXYLIC_ACID
            is_acidic = True
        elif not is_nitrosamine and _BASE_SMARTS and mol.HasSubstructMatch(_BASE_SMARTS):
            pka = _PKA_ALIPHATIC_AMINE
            is_basic = True
        elif (
            not is_nitrosamine
            and _AROMATIC_AMINE_SMARTS
            and mol.HasSubstructMatch(_AROMATIC_AMINE_SMARTS)
        ):
            # 아릴아민 — pH 7.4에서 대부분 중성이라 LogD ≈ LogP
            pka = _PKA_AROMATIC_AMINE
            is_basic = True

        if is_acidic:
            logd = logp_baseline - math.log10(1 + 10 ** (target_ph - pka))
        elif is_basic:
            logd = logp_baseline - math.log10(1 + 10 ** (pka - target_ph))
        else:
            logd = logp_baseline

        return {
            "MW": round(mw, 2),
            "LogP_Baseline": round(logp_baseline, 2),
            "LogD_pH7.4": round(logd, 2),
            "Estimated_pKa": round(pka, 2),
            "Boiling_Point_C": round(boiling_point, 1),
            "TPSA": round(tpsa, 2),
        }

    @staticmethod
    def evaluate_thermodynamic_viability(props: dict[str, float] | None) -> bool:
        if not props:
            return False
        if props.get("Boiling_Point_C", 0) < -100.0:
            return False
        if props.get("MW", 0) > 1000.0:
            return False
        return True


class DualLayerReactionParser:
    """Main vs side parallel SMIRKS parsing — delegates to process_reaction_rules."""

    def parse_process_matrix(self, reactant_smiles_list: list[str]) -> dict[str, Any]:
        from p31_process_reaction_rules import predict_process_reactions
        from p91_prime_integrated_report import format_process_synthesis_analysis

        mols_ok = [
            s
            for s in reactant_smiles_list
            if s and Chem.MolFromSmiles(s.strip())
        ]
        if len(mols_ok) < 2:
            target = mols_ok[0] if mols_ok else None
            return {"target": target, "process_impurities": []}

        raw = predict_process_reactions(reactant_smiles_list)
        syn = format_process_synthesis_analysis(raw)

        impurities = [
            {
                "smiles": row["Side_Product_SMILES"],
                "origin": f"공정 부반응식 ({syn.get('Triggered_Mechanism') or 'Side'})",
                "alert": row.get("Alert_Description", ""),
            }
            for row in syn.get("Process_Impurities") or []
        ]

        return {
            "target": syn.get("Main_Product_SMILES") or None,
            "process_impurities": impurities,
            "triggered_mechanism": syn.get("Triggered_Mechanism", ""),
        }


class DegradationTreeEngine:
    """
    Recursive DFS degradation tree — ICH Q3A/B in-vitro shelf-life degradants only.

    - `self.rules`: 가수분해 + 비라디칼 자동산화 SMIRKS (`tree_engine_degradation_rules`)
    - Child 여과: `ShelfLifeChildFilter` (라디칼 · 생체 대사 배제)
    """

    def __init__(self) -> None:
        from p41_degradation_smirks_rules import tree_engine_degradation_rules

        self.rules = {
            name: info["smirks"] for name, info in tree_engine_degradation_rules().items()
        }

    def _canonical(self, smiles: str) -> str | None:
        canon, _ = safe_canonical_smiles(smiles)
        return canon

    def _passes_shelf_life_child_filter(
        self,
        parent_mol: Chem.Mol | None,
        child_mol: Chem.Mol | None,
        *,
        parent_smiles: str = "",
        child_smiles: str = "",
    ) -> bool:
        """라디칼 · `[OX2H]-[c;r6]` 방향족 수산화 정밀 배제."""
        try:
            return ShelfLifeChildFilter.evaluate_mols(
                parent_mol,
                child_mol,
                parent_smiles=parent_smiles,
                child_smiles=child_smiles,
            ).eligible
        except Exception:
            return True

    @staticmethod
    def _is_trackable_fragment(smiles: str, *, parent: str) -> bool:
        if not smiles or smiles == parent:
            return False
        if not is_spin_stabilized_neutral_smiles(smiles):
            return False
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        if Descriptors.MolWt(mol) < _MIN_DEGRADANT_MW:
            return False
        if mol.GetNumHeavyAtoms() < 2:
            return False
        if _NITROSO_FRAGMENT_SMARTS and mol.HasSubstructMatch(_NITROSO_FRAGMENT_SMARTS):
            if not (_NITROSO_SMARTS and mol.HasSubstructMatch(_NITROSO_SMARTS)):
                return False
        return True

    def _run_rule(self, mol: Chem.Mol, smirks: str, *, parent: str) -> list[str]:
        smirks = (smirks or "").replace(" ", "")
        try:
            rxn = AllChem.ReactionFromSmarts(smirks)
        except Exception:
            return []
        if rxn is None:
            return []

        out: list[str] = []
        seen: set[str] = set()
        try:
            product_sets = rxn.RunReactants((mol,))
        except Exception:
            return []

        for prod_tuple in product_sets:
            for prod_mol in prod_tuple:
                try:
                    if not is_spin_stabilized_neutral_mol(prod_mol):
                        continue
                    Chem.SanitizeMol(prod_mol)
                    # 부모에 없던 [OX2H]-[c;r6] → 대사산물로 제외
                    if should_cull_metabolic_aromatic_hydroxylation(mol, prod_mol):
                        continue
                    smi = Chem.MolToSmiles(prod_mol, canonical=True, isomericSmiles=True)
                    if not smi or smi in seen:
                        continue
                    if not self._is_trackable_fragment(smi, parent=parent):
                        continue
                    seen.add(smi)
                    out.append(smi)
                except Exception:
                    continue
        return out

    def generate_tree(
        self,
        smiles: str,
        current_depth: int = 1,
        max_depth: int = 3,
        visited_set: set[str] | None = None,
        *,
        formulation_ph: float = 7.4,
    ) -> dict[str, Any] | None:
        if visited_set is None:
            visited_set = set()

        canon = self._canonical(smiles)
        if not canon or canon in visited_set or current_depth > max_depth:
            return None

        local_visited = visited_set.copy()
        local_visited.add(canon)

        props = HybridPropertyEngine.calculate_properties(canon, target_ph=formulation_ph)
        if not HybridPropertyEngine.evaluate_thermodynamic_viability(props):
            return None

        node: dict[str, Any] = {
            "smiles": canon,
            "depth": current_depth,
            "properties": props,
            "mechanism_origin": "Parent Core" if current_depth == 1 else "Degradation Split",
            "children": [],
        }

        mol = Chem.MolFromSmiles(canon)
        if mol is None or current_depth >= max_depth:
            return node
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return node

        child_seen: set[str] = set()
        for rule_name, smirks in self.rules.items():
            for child_smiles in self._run_rule(mol, smirks, parent=canon):
                if child_smiles == canon or child_smiles in child_seen:
                    continue

                child_canon = self._canonical(child_smiles) or child_smiles
                child_mol = Chem.MolFromSmiles(child_canon)
                if child_mol is None:
                    continue
                try:
                    Chem.SanitizeMol(child_mol)
                except Exception:
                    continue

                # canonical SMILES 도출 직후 · 재귀 호출 직전 — 정밀 여과
                if not self._passes_shelf_life_child_filter(
                    mol,
                    child_mol,
                    parent_smiles=canon,
                    child_smiles=child_canon,
                ):
                    continue

                child_seen.add(child_canon)
                child_tree = self.generate_tree(
                    child_canon,
                    current_depth + 1,
                    max_depth,
                    local_visited,
                    formulation_ph=formulation_ph,
                )
                if child_tree:
                    child_tree["mechanism_origin"] = rule_name
                    node["children"].append(child_tree)

        return node


def _ames_display_label(ames_result: str) -> str:
    val = (ames_result or "unknown").strip().lower()
    if val == "positive":
        return "🔥 Positive (실측)"
    if val == "negative":
        return "🟢 Negative (실측)"
    return "⚪ Unknown"


# 2026-09-03 — 실측 DB를 못 읽었을 때 붙는 출처 라벨.
#   전에는 fallback 항목도 `source: "Hansen"` 이라 보고서에서 실측 풀과
#   구분되지 않았습니다.
HANSEN_FALLBACK_SOURCE = (
    "Hansen (built-in 7-compound fallback — empirical DB unavailable)"
)


def _load_read_across_pool() -> tuple[list[dict[str, str]], bool, str]:
    """read-across 풀을 읽습니다.

    반환: (풀, degraded 여부, 사유)

    2026-09-03 수정 — 전에는 리스트만 돌려주고 실패를
    `except Exception: return copy.deepcopy(HANSEN_MUTAGENICITY_DB)` 로
    삼켰습니다. `_ensure_db()`가 FileNotFoundError를 던지면 풀이 조용히
    **7개 화합물**로 줄어드는데, docstring과 보고서 라벨은 계속
    "Hansen/ISSCAN/Tox21 empirical pool"을 주장하고 `find_analogues`는
    유사체 0건을 정상 결과처럼 반환했습니다. 사용자는 read-across가
    수행되지 않았다는 사실을 알 수 없었습니다.
    이제 degraded 플래그와 사유를 함께 올려 보고서에 노출시킵니다.
    """
    try:
        from p60_m7_empirical_db import _ensure_db, get_default_db_connection

        _ensure_db()
        conn = get_default_db_connection()
        rows = conn.execute(
            """
            SELECT compound_name, canonical_smiles, ames_result, source
            FROM m7_empirical
            """
        ).fetchall()
        if not rows:
            raise ValueError("empty empirical db")

        by_smi: dict[str, dict[str, str]] = {}
        rank = {"positive": 3, "negative": 2, "unknown": 1}
        for row in rows:
            smi = str(row["canonical_smiles"])
            ames = str(row["ames_result"] or "unknown").lower()
            entry = by_smi.get(smi)
            if entry is None or rank.get(ames, 0) > rank.get(entry["ames_result"], 0):
                by_smi[smi] = {
                    "name": str(row["compound_name"] or smi),
                    "smiles": smi,
                    "ames_result": ames,
                    "source": str(row["source"] or "Hansen"),
                }
        return list(by_smi.values()), False, ""
    except Exception as exc:
        import logging

        reason = f"{type(exc).__name__}: {exc}"
        logging.getLogger(__name__).warning(
            "실측 read-across DB를 읽지 못했습니다 (%s) — 내장 7종 fallback으로 "
            "축소됩니다. 유사체 0건은 «해당 없음»이 아니라 «검색 풀이 축소됨»을 "
            "뜻할 수 있습니다.",
            reason,
        )
        pool = copy.deepcopy(HANSEN_MUTAGENICITY_DB)
        for item in pool:
            item["source"] = HANSEN_FALLBACK_SOURCE
        return pool, True, reason


def _merge_read_across_extras(pool: list[dict[str, str]]) -> list[dict[str, str]]:
    extras = [
        {
            "name": "N-Nitrosodiethylamine (NDEA)",
            "smiles": "CCN(CC)N=O",
            "ames_result": "positive",
            "source": "Hansen",
        },
    ]
    seen = {item["smiles"] for item in pool}
    out = list(pool)
    for row in extras:
        if row["smiles"] not in seen:
            out.append(row)
            seen.add(row["smiles"])
    return out


def read_across_pool_status() -> dict[str, Any]:
    """현재 read-across 풀의 상태 — 보고서에 그대로 실을 수 있는 형태.

    2026-09-03 추가. `find_analogues`가 유사체 0건을 돌려줄 때 그것이
    «구조적으로 유사체가 없다»인지 «검색 풀이 7종으로 축소됐다»인지
    구분할 수 있어야 합니다.
    """
    pool, degraded, reason = _load_read_across_pool()
    pool = _merge_read_across_extras(pool)
    return {
        "pool_size": len(pool),
        "pool_degraded": degraded,
        "pool_degraded_reason": reason,
        "pool_source": (
            HANSEN_FALLBACK_SOURCE
            if degraded
            else "Hansen/ISSCAN/Tox21 empirical pool (SQLite)"
        ),
    }


class ReadAcrossEngine:
    """Morgan Tanimoto read-across against Hansen/ISSCAN/Tox21 empirical pool.

    ⚠️ 실측 SQLite 인덱스를 읽지 못하면 내장 7종 fallback으로 축소됩니다.
       그 경우 결과 각 행과 `find_analogues_with_pool_info()`의
       `pool_degraded` 플래그로 그 사실이 드러납니다.
    """

    @staticmethod
    def find_analogues(
        target_smiles: str,
        top_n: int = 3,
        *,
        min_similarity: float = READ_ACROSS_MIN_SIMILARITY,
    ) -> list[dict[str, Any]]:
        """유사체 목록. 반환 형식(리스트)은 하위 호환을 위해 유지합니다.

        2026-09-03 — 각 행에 `pool_size`·`pool_degraded`를 실어
        보고서가 «풀이 축소된 상태에서 나온 결과»임을 알 수 있게 했습니다.
        0건이라 실을 곳이 없는 경우를 위해
        `find_analogues_with_pool_info()`와 `read_across_pool_status()`를
        따로 두었습니다.
        """
        return ReadAcrossEngine.find_analogues_with_pool_info(
            target_smiles, top_n, min_similarity=min_similarity
        )["analogues"]

    @staticmethod
    def find_analogues_with_pool_info(
        target_smiles: str,
        top_n: int = 3,
        *,
        min_similarity: float = READ_ACROSS_MIN_SIMILARITY,
    ) -> dict[str, Any]:
        """유사체 + 검색 풀 상태를 함께 돌려줍니다.

        반환: {"analogues", "pool_size", "pool_degraded",
               "pool_degraded_reason", "pool_source"}
        """
        raw_pool, degraded, reason = _load_read_across_pool()
        pool = _merge_read_across_extras(raw_pool)
        pool_source = (
            HANSEN_FALLBACK_SOURCE
            if degraded
            else "Hansen/ISSCAN/Tox21 empirical pool (SQLite)"
        )
        envelope: dict[str, Any] = {
            "analogues": [],
            "pool_size": len(pool),
            "pool_degraded": degraded,
            "pool_degraded_reason": reason,
            "pool_source": pool_source,
        }

        target_mol = Chem.MolFromSmiles((target_smiles or "").strip())
        if target_mol is None:
            return envelope

        target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, 2, nBits=2048)
        results: list[dict[str, Any]] = []
        min_sim = max(0.0, min(1.0, float(min_similarity)))

        for item in pool:
            db_mol = Chem.MolFromSmiles(item["smiles"])
            if db_mol is None:
                continue
            db_fp = AllChem.GetMorganFingerprintAsBitVect(db_mol, 2, nBits=2048)
            similarity = DataStructs.TanimotoSimilarity(target_fp, db_fp)
            if similarity < min_sim:
                continue
            results.append(
                {
                    "name": item["name"],
                    "smiles": item["smiles"],
                    "ames_result": _ames_display_label(item["ames_result"]),
                    "ames_raw": item["ames_result"],
                    "source": item.get("source", "Hansen"),
                    "similarity": round(similarity, 3),
                    # 이 결과가 어떤 풀에서 나왔는지 행마다 남깁니다
                    "pool_size": len(pool),
                    "pool_degraded": degraded,
                }
            )

        results.sort(key=lambda x: x["similarity"], reverse=True)
        envelope["analogues"] = results[:top_n]
        return envelope


def _export_read_across_name(name: str, similarity: float) -> str:
    n = (name or "").strip()
    if similarity >= 0.999:
        for suffix in (" (NDEA)", " (NDMA)", " (MNU)"):
            if n.endswith(suffix):
                return n[: -len(suffix)]
    return n


def format_read_across_export(
    hits: list[dict[str, Any]],
    *,
    node_display_smiles: str | None = None,
    node_canonical: str | None = None,
    include_structure_b64: bool = True,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    node_canon = (node_canonical or "").strip()
    display_smi = (node_display_smiles or "").strip()

    for hit in hits:
        hit_smi = str(hit.get("smiles") or "")
        hit_canon, _ = safe_canonical_smiles(hit_smi)
        smiles_out = hit_smi
        if node_canon and hit_canon == node_canon and display_smi:
            smiles_out = display_smi

        sim = float(hit.get("similarity", 0))
        row: dict[str, Any] = {
            "name": _export_read_across_name(str(hit.get("name") or ""), sim),
            "smiles": smiles_out,
            "ames_result": hit["ames_result"],
            "similarity": hit["similarity"],
        }
        # 2026-09-03: 축소된 풀에서 나온 결과임을 보고서까지 전달합니다
        if "pool_size" in hit:
            row["pool_size"] = hit["pool_size"]
        if hit.get("pool_degraded"):
            row["pool_degraded"] = True
            row["source"] = HANSEN_FALLBACK_SOURCE
        if include_structure_b64:
            b64 = smiles_to_base64_img(smiles_out, size=(200, 200))
            if b64:
                row["structure_image_base64"] = b64
        rows.append(row)
    return rows


def format_properties_export(props: dict[str, float] | None) -> dict[str, float]:
    props = props or {}

    def _r(key: str, ndigits: int = 2) -> float:
        try:
            return round(float(props.get(key, 0.0)), ndigits)
        except (TypeError, ValueError):
            return 0.0

    return {
        "MW": _r("MW", 2),
        "LogD_pH7.4": _r("LogD_pH7.4", 2),
        "Estimated_pKa": _r("Estimated_pKa", 2),
        "Boiling_Point_C": _r("Boiling_Point_C", 1),
    }


def enrich_tree_for_api(
    tree: dict[str, Any] | None,
    *,
    formulation_ph: float = 7.4,
    read_across_top_n: int = 1,
    min_similarity: float = READ_ACROSS_MIN_SIMILARITY,
    root_display_smiles: str | None = None,
    include_structure_b64: bool = True,
) -> dict[str, Any] | None:
    if not tree:
        return None

    canon_smiles = str(tree.get("smiles") or "")
    node_canon, _ = safe_canonical_smiles(canon_smiles)
    display_smiles = canon_smiles
    if tree.get("depth") == 1 and root_display_smiles:
        display_smiles = root_display_smiles

    # 2026-09-03: 유사체가 0건일 때도 «풀이 축소됐는지»를 알 수 있도록
    #   envelope 형태로 받아 노드에 플래그를 함께 싣습니다.
    ra = ReadAcrossEngine.find_analogues_with_pool_info(
        canon_smiles or display_smiles,
        top_n=read_across_top_n,
        min_similarity=min_similarity,
    )
    analogues = ra["analogues"]

    node: dict[str, Any] = {
        "smiles": display_smiles,
        "depth": tree.get("depth", 1),
        "mechanism_origin": tree.get("mechanism_origin", "Parent Core"),
        "properties": format_properties_export(tree.get("properties")),
        "read_across_analogues": format_read_across_export(
            analogues,
            node_display_smiles=display_smiles,
            node_canonical=node_canon or canon_smiles,
            include_structure_b64=include_structure_b64,
        ),
        "children": [
            child
            for child in (
                enrich_tree_for_api(
                    c,
                    formulation_ph=formulation_ph,
                    read_across_top_n=read_across_top_n,
                    min_similarity=min_similarity,
                    include_structure_b64=include_structure_b64,
                )
                for c in (tree.get("children") or [])
            )
            if child
        ],
        # read-across 검색 풀 상태 — 유사체 0건의 의미를 구분할 수 있게 합니다
        "read_across_pool_size": ra["pool_size"],
        "read_across_pool_degraded": ra["pool_degraded"],
        "read_across_pool_source": ra["pool_source"],
    }
    if ra["pool_degraded"] and ra.get("pool_degraded_reason"):
        node["read_across_pool_degraded_reason"] = ra["pool_degraded_reason"]
    if include_structure_b64:
        node_b64 = smiles_to_base64_img(display_smiles, size=(200, 200))
        if node_b64:
            node["structure_image_base64"] = node_b64
    return node


def _display_parent_smiles(user_input: str, canonical: str) -> str:
    """Preserve user-entered SMILES notation when chemically equivalent."""
    raw = (user_input or "").strip()
    if not raw:
        return canonical
    raw_canon, err = safe_canonical_smiles(raw)
    if not err and raw_canon == canonical:
        return raw
    return canonical


def build_degradation_api_response(
    parent_smiles: str,
    *,
    formulation_ph: float = 7.4,
    max_depth: int = 3,
    read_across_top_n: int = 1,
    min_similarity: float = READ_ACROSS_MIN_SIMILARITY,
    include_structure_b64: bool = True,
) -> dict[str, Any]:
    """
    Frontend / regulatory API payload:
    status, parent_molecule, formulation_ph, degradation_tree (+ per-node read-across).
    """
    parent = (parent_smiles or "").strip()
    canon, err = safe_canonical_smiles(parent)
    if err or not canon:
        return {
            "status": "error",
            "parent_molecule": parent,
            "formulation_ph": formulation_ph,
            "message": err or "Invalid SMILES",
            "degradation_tree": None,
        }

    raw_tree = DegradationTreeEngine().generate_tree(
        canon,
        max_depth=max_depth,
        formulation_ph=formulation_ph,
    )
    if not raw_tree:
        return {
            "status": "error",
            "parent_molecule": canon,
            "formulation_ph": formulation_ph,
            "message": "Thermodynamic filter rejected parent or invalid structure",
            "degradation_tree": None,
        }

    display_parent = _display_parent_smiles(parent, canon)
    tree_api = enrich_tree_for_api(
        raw_tree,
        formulation_ph=formulation_ph,
        read_across_top_n=read_across_top_n,
        min_similarity=min_similarity,
        root_display_smiles=display_parent,
        include_structure_b64=include_structure_b64,
    )

    return {
        "status": "success",
        "parent_molecule": display_parent,
        "formulation_ph": formulation_ph,
        "degradation_tree": tree_api,
    }


def degradation_api_response_to_json(payload: dict[str, Any], *, indent: int = 2) -> str:
    import json

    return json.dumps(payload, indent=indent, ensure_ascii=False)


def run_prime_backend_pipeline(
    reactant_smiles_list: list[str],
    degradation_parent_smiles: str | None = None,
    *,
    max_depth: int = 3,
    read_across_top_n: int = 3,
) -> dict[str, Any]:
    """Unified backend demo / API payload."""
    parser = DualLayerReactionParser()
    process = parser.parse_process_matrix(reactant_smiles_list)

    parent = (degradation_parent_smiles or process.get("target") or "").strip()
    tree = None
    if parent:
        tree = DegradationTreeEngine().generate_tree(parent, max_depth=max_depth)

    analogues: list[dict[str, Any]] = []
    ra_target = parent or (reactant_smiles_list[0] if reactant_smiles_list else "")
    if ra_target:
        analogues = ReadAcrossEngine.find_analogues(ra_target, top_n=read_across_top_n)

    return {
        "Process_Dual_Layer": process,
        "Degradation_Tree_With_Properties": tree,
        "Read_Across_Analogues": analogues,
    }


if __name__ == "__main__":
    import json

    print("=" * 80)
    print(" [PRIME-DS Backend Core Suite] 초고도화 엔진 가동 세션")
    print("=" * 80)

    ndea_api = build_degradation_api_response("CCN(N=O)CC", formulation_ph=7.4, max_depth=2)
    print("\n[API] NDEA enriched degradation tree:")
    print(degradation_api_response_to_json(ndea_api))

    reactants = ["CC(=O)O", "CCNCC"]
    result = run_prime_backend_pipeline(
        reactants,
        degradation_parent_smiles="CCOC(=O)CCSCCc1ccccc1",
        max_depth=2,
    )

    proc = result["Process_Dual_Layer"]
    print("\n[1] 공정 Dual-Layer 파싱 결과:")
    print(f"  - 생성 타겟 의약품(Target Core): {proc.get('target')}")
    print(f"  - 메커니즘: {proc.get('triggered_mechanism')}")
    print("  - 검출된 경쟁 공정 불순물(Side Matrix):")
    for imp in proc.get("process_impurities") or []:
        print(f"    ▶ 구조: {imp['smiles']} | 출처: {imp['origin']}")

    tree = result["Degradation_Tree_With_Properties"]
    print("\n[2] 재귀적 다단계 분해 (Thermodynamic Viable Tree):")
    if tree:
        print(f"  - Parent: {tree['smiles']} | BP={tree['properties']['Boiling_Point_C']}°C")
        for child in tree.get("children") or []:
            print(
                f"    └ Gen-{child['depth']} {child['smiles']} "
                f"({child['mechanism_origin']}) LogD={child['properties']['LogD_pH7.4']}"
            )

    print("\n[3] Read-Across (Hansen/ISSCAN Tanimoto Top-3):")
    for hit in result["Read_Across_Analogues"]:
        print(
            f"  - {hit['name']} | Tanimoto={hit['similarity']} | "
            f"{hit['ames_result']} | {hit['source']}"
        )
