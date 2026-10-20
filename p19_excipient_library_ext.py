"""부형제 라이브러리 확장 — 구조 보강 + 신규 등재.

두 가지를 합니다.

① **구조 보강** (`STRUCTURE_PATCH`)
   기존 76종 중 25종은 `smiles`도 `reactive_probes`도 없어 반응 예측에
   전혀 참여하지 못했습니다. 골라도 항상 «이상 없음»이 나옵니다.
   메타중아황산나트륨·젤라틴·레시틴처럼 실제 반응성이 큰 것들이
   포함되어 있어 위음성입니다. 여기에 구조와 반응 조각을 채웁니다.

② **신규 등재** (`ADDITIONAL_EXCIPIENTS`)
   경구고형제에 흔히 쓰이나 목록에 없던 부형제를 추가합니다.

--------------------------------------------------------------------------
`reactive_probes`가 핵심입니다
--------------------------------------------------------------------------
고분자·무기물은 전체 구조를 SMILES로 쓸 수 없습니다. 대신 **실제로
반응하는 조각**을 적습니다. 유당이 좋은 예입니다 — 고리 구조로는
Maillard가 안 걸리지만, 용액에서 열린 알데하이드 조각을 프로브로 두면
정확히 걸립니다.

    "Lactose_Monohydrate": {
        "smiles": "OCC1OC(...)C1O",        # 전체 구조 (닫힌 고리)
        "reactive_probes": ["O=CC(O)CO"],  # 반응하는 조각 (열린 알데하이드)
    }

고분자에는 **반복 단위 또는 유리되는 분해산물**을 프로브로 씁니다.
예: 크로스포비돈 → 잔류 과산화물 `OO`, CAP → 아세트산 `CC(=O)O`.

--------------------------------------------------------------------------
출처
--------------------------------------------------------------------------
구조는 PubChem canonical SMILES, 명칭은 USP-NF / Ph.Eur. 표제어를
따랐습니다. 배합금기 문구는 Handbook of Pharmaceutical Excipients의
Incompatibilities 항목을 요약한 것입니다.

⚠️ 고분자의 `smiles`는 **대표 단량체 또는 대표 조각**입니다. 실제 고분자
   구조가 아니며, 몰비 계산에 쓰면 안 됩니다 (화면에서 «정의 없음»으로
   표시됩니다).

--------------------------------------------------------------------------
2026-09-03 개정
--------------------------------------------------------------------------
① `mw_kind` 판정표(`MW_KIND_BY_KEY`)를 추가했습니다. 단량체 대용 SMILES도
   RDKit에서는 멀쩡한 분자량을 반환하기 때문에, 위 «몰비에 쓰면 안 된다»는
   주의가 **코드로는 지켜지지 않고 있었습니다**. p17의 `if not exc_mw`
   가드가 한 번도 발동하지 않았습니다. 이제 항목마다 명시적으로 표시합니다.

② 구조 오류 세 건을 고쳤습니다.
     Sodium_Metabisulfite          Na2S2O4(dithionite) → Na2S2O5  (MW 174.11→190.11)
     Carboxymethylcellulose_Sodium 순전하 −1 → [Na+] 추가
     Eudragit_RS                   순전하 +1 → [Cl-] 추가

③ Isopropyl_Alcohol의 반응 프로브가 아세트알데하이드(에탄올 유래)였습니다.
   IPA는 2급 알코올이라 산화하면 아세톤이 됩니다 → `CC(C)=O`.

④ `validate_extension_smiles()` — 전수 파싱·전하·mw_kind 검증.

2026-08-25 작성.
"""

from __future__ import annotations

from typing import Any

__all__ = [
    "STRUCTURE_PATCH",
    "PROBE_PATCH",
    "ADDITIONAL_EXCIPIENTS",
    "MW_KIND_BY_KEY",
    "apply_extensions",
    "extension_stats",
    "validate_extension_smiles",
]


# ===========================================================================
# 반응 조각 보충 — 이미 프로브가 있지만 빠진 것을 더합니다
#
# 같은 계열인데 프로브가 서로 다른 경우가 있습니다. 확립된 과산화물
# 공급원 계열을 전수 확인한 결과 **PEG 4000 하나만** 빠져 있었습니다.
#
#   크로스포비돈 · 포비돈 K30/K90        ['C(=O)O', 'OO']  ✅
#   폴리소르베이트 20/60/80 · 폴록사머 407 ['C(=O)O', 'OO']  ✅
#   PEG 4000                           ['C(=O)O', 'O=CO'] ❌ 과산화물 없음
#
# 마크로골은 자동산화로 과산화물이 쌓이는 대표 부형제인데 이것만
# 빠져 있어, 3급 아민만 가진 원료(아미오다론 등)와 조합하면
# «상호작용 없음»이 나왔습니다.
#
# ⚠️ 셀룰로오스 에터(HPMC·HPC·EC)는 넣지 않았습니다.
#    미량 과산화물이 보고되기는 하지만 포비돈·마크로골·폴리소르베이트만큼
#    확립된 경로가 아니고, 넣으면 거의 모든 고분자 부형제가 3급 아민
#    원료와 걸려 변별력이 떨어집니다.
# ===========================================================================

PROBE_PATCH: dict[str, list[str]] = {
    # 마크로골 — 자동산화 과산화물 + 아세트알데하이드
    # (포름산 'C(=O)O'는 이미 등록되어 있습니다)
    "Polyethylene_Glycol_4000": ["OO", "O=CC"],
}


# ===========================================================================
# ① 구조 보강 — 기존 항목에 smiles / reactive_probes 를 채웁니다
# ===========================================================================

STRUCTURE_PATCH: dict[str, dict[str, Any]] = {
    # ── 단백질 · 다당류 ──────────────────────────────────────────────────
    "Gelatin": {
        # 대표 반복 단위 (Gly-Pro-Hyp) — 유리 아미노기가 반응 자리
        "smiles": "NCC(=O)N1CCCC1C(=O)NCC(=O)O",
        "reactive_probes": ["NCC(=O)O", "NCCCCC(N)C(=O)O"],
    },
    "Xanthan_Gum": {
        "smiles": "OCC1OC(O)C(O)C(O)C1O",
        "reactive_probes": ["OC(=O)C1OC(O)C(O)C(O)C1O"],
    },
    "Carboxymethylcellulose_Sodium": {
        # 2026-09-03 수정: 카복실레이트 음이온만 있고 대이온 Na⁺가 빠져
        #   순전하 −1 (C8H13O8⁻)이었습니다. 나트륨염이므로 [Na+]를 더해
        #   중성 C8H13NaO8 (MW 260.17)로 맞춥니다.
        "smiles": "[Na+].OCC1OC(O)C(OCC(=O)[O-])C(O)C1O",
        "reactive_probes": ["OCC(=O)O"],
    },
    "Alginic_Acid": {
        "smiles": "OC1C(O)C(OC(C1O)C(=O)O)O",
        "reactive_probes": ["OC1C(O)C(O)OC(C(=O)O)C1O"],
    },
    "Acacia_Gum": {
        "smiles": "OCC1OC(O)C(O)C(O)C1O",
        "reactive_probes": ["O=CC(O)C(O)C(O)CO"],
    },

    # ── 무기염 · 완충제 ──────────────────────────────────────────────────
    # 반응 조각이 없는 것이 정상입니다. 구조만 채워 몰비를 계산할 수 있게 합니다.
    "Sodium_Chloride": {"smiles": "[Na+].[Cl-]", "reactive_probes": []},
    "Calcium_Phosphate_Tribasic": {
        "smiles": "[Ca+2].[Ca+2].[Ca+2].[O-]P([O-])([O-])=O.[O-]P([O-])([O-])=O",
        "reactive_probes": [],
    },
    "Titanium_Dioxide": {"smiles": "O=[Ti]=O", "reactive_probes": []},
    "Sodium_Phosphate_Dibasic": {
        "smiles": "[Na+].[Na+].OP([O-])([O-])=O",
        "reactive_probes": [],
    },
    "Potassium_Phosphate": {
        "smiles": "[K+].[K+].OP([O-])([O-])=O",
        "reactive_probes": [],
    },
    "Sodium_Hydroxide": {"smiles": "[Na+].[OH-]", "reactive_probes": []},
    "Hydrochloric_Acid": {"smiles": "Cl", "reactive_probes": []},

    # ── 지질 · 유지 ──────────────────────────────────────────────────────
    "Glyceryl_Monostearate": {
        "smiles": "CCCCCCCCCCCCCCCCCC(=O)OCC(O)CO",
        # 에스터 → 에스터 교환 · 아미놀리시스, 가수분해 시 스테아르산
        "reactive_probes": ["CCCCCCCCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCCCCCC(=O)O"],
    },
    "Hydrogenated_Castor_Oil": {
        "smiles": "CCCCCCC(O)CCCCCCCCCCC(=O)OCC(O)CO",
        "reactive_probes": ["CCCCCCC(O)CCCCCCCCCCC(=O)O"],
    },
    "Light_Mineral_Oil": {
        "smiles": "CCCCCCCCCCCCCCCC",
        "reactive_probes": [],
    },
    "Lecithin": {
        # 포스파티딜콜린 — 에스터 결합 + 4급 암모늄
        "smiles": "CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC",
        "reactive_probes": ["CCCCCCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCCCC(=O)O"],
    },
    "Benzyl_Benzoate": {
        "smiles": "O=C(OCc1ccccc1)c1ccccc1",
        "reactive_probes": ["O=C(OC)c1ccccc1", "O=C(O)c1ccccc1"],
    },

    # ── 방부제 · 항산화제 ────────────────────────────────────────────────
    "Sodium_Benzoate": {
        "smiles": "[Na+].[O-]C(=O)c1ccccc1",
        "reactive_probes": ["OC(=O)c1ccccc1"],
    },
    "Potassium_Sorbate": {
        # α,β-불포화 → Michael 수용체
        "smiles": "[K+].C/C=C/C=C/C(=O)[O-]",
        "reactive_probes": ["C/C=C/C=C/C(=O)O"],
    },
    "Butylated_Hydroxyanisole_BHA": {
        "smiles": "COc1ccc(O)c(C(C)(C)C)c1",
        "reactive_probes": ["COc1ccc(O)c(C(C)(C)C)c1"],
    },
    "Sodium_Metabisulfite": {
        # ⚠️ 아황산염 — 알데하이드·케톤 API와 중아황산 부가물 형성
        #
        # 2026-09-03 수정: 기존 SMILES `[O-]S(=O)S([O-])=O` 는 **S–S 결합**이라
        #   Na2S2O4 (아이티온산나트륨 / dithionite, 실측 Na2O4S2 · MW 174.11)
        #   였습니다. 메타중아황산나트륨은 S–O–S 로 이어진 Na2S2O5 (MW 190.11)
        #   입니다. 16 g/mol 틀린 값이 몰비 계산에 들어가고 있었습니다.
        #   교체 후 실측 Na2O5S2 · MW 190.11 확인.
        "smiles": "[Na+].[Na+].[O-]S(=O)OS(=O)[O-]",
        "reactive_probes": ["OS(=O)O", "[O-]S(=O)O"],
    },

    # ── 장용 코팅 ────────────────────────────────────────────────────────
    "Hypromellose_Phthalate": {
        "smiles": "OC(=O)c1ccccc1C(=O)OCC1OC(O)C(O)C(O)C1O",
        "reactive_probes": ["OC(=O)c1ccccc1C(=O)O", "OC(=O)c1ccccc1C(=O)OC"],
    },
    "Cellulose_Acetate_Phthalate": {
        # 가수분해로 아세트산 · 프탈산 유리
        "smiles": "CC(=O)OCC1OC(O)C(OC(=O)c2ccccc2C(=O)O)C(O)C1O",
        "reactive_probes": ["CC(=O)O", "OC(=O)c1ccccc1C(=O)O", "CC(=O)OC"],
    },

    # ── 감미료 ───────────────────────────────────────────────────────────
    "Saccharin_Sodium": {
        "smiles": "[Na+].O=C1[N-]S(=O)(=O)c2ccccc21",
        "reactive_probes": [],
    },

    # ── 프로브는 있는데 구조가 비었거나 무효한 것들 ──────────────────
    #    반응 예측은 되지만 분자량을 못 구해 몰비가 «정의 없음»으로
    #    나오고 있었습니다.
    "Povidone_K90": {
        # 대표 반복 단위(N-비닐피롤리돈). K30과 같은 값입니다.
        "smiles": "CCN1CCCC1=O",
        "reactive_probes": ["C(=O)O", "OO"],
    },
    "Dextrose": {
        "smiles": "OCC1OC(O)C(O)C(O)C1O",
        "reactive_probes": ["O=CC(O)CO"],
    },
    "Isopropyl_Alcohol": {
        # 2026-09-03 수정: 프로브가 아세트알데하이드 `CC=O` 였습니다.
        #   그건 **에탄올**의 산화 산물입니다. IPA는 2급 알코올이라 산화하면
        #   아세톤 `CC(C)=O` 가 됩니다. 알데하이드는 Maillard 이민을 만들고
        #   케톤은 만들지 않으므로, 기존 값은 **위양성** 방향이었습니다.
        "smiles": "CC(C)O",
        "reactive_probes": ["CC(C)=O"],
    },
    "Polysorbate_60": {
        "smiles": "CCCCCCCCCCCCCCCCCC(=O)OCCOCCO",
        "reactive_probes": ["C(=O)O", "OO"],
    },
    "Magnesium_Stearate": {
        # ⚠️ 원본 SMILES가 문법 오류였습니다 —
        #    두 번째 스테아레이트 사슬 전체가 대괄호로 묶여 있었습니다.
        #        CCCCCCCCCCCCCCCCCC(=O)[O-].[CCCCCCCCCCCCCCCCCC(=O)[O-]].[Mg+2]
        #                                    ↑ 여기
        #    RDKit이 읽지 못해 몰비가 나오지 않았습니다.
        "smiles": "CCCCCCCCCCCCCCCCCC(=O)[O-].CCCCCCCCCCCCCCCCCC(=O)[O-].[Mg+2]",
        "reactive_probes": ["CCCCCCCCCCCCCCCCCC(=O)O"],
    },
    "Peppermint_Oil": {
        # 대표 성분 멘톨 + 멘톤(케톤 → 중아황산 부가 가능)
        "smiles": "CC(C)C1CCC(C)CC1O",
        "reactive_probes": ["CC(C)C1CCC(C)CC1=O"],
    },
}


# ===========================================================================
# ② 신규 등재
# ===========================================================================

ADDITIONAL_EXCIPIENTS: dict[str, dict[str, Any]] = {
    # ── Binder / Diluent ─────────────────────────────────────────────────
    "Maltodextrin": {
        "smiles": "OCC1OC(O)C(O)C(O)C1O",
        "category": "Binder/Diluent",
        "label_ko": "말토덱스트린 (Maltodextrin)",
        "reactive_probes": ["O=CC(O)C(O)C(O)C(O)CO"],
        "alerts": "🔥 환원 말단 — 1급 아민 API Maillard",
    },
    "Maltose": {
        "smiles": "OCC1OC(OC2C(O)C(O)C(O)OC2CO)C(O)C(O)C1O",
        "category": "Binder/Diluent",
        "label_ko": "맥아당 (Maltose)",
        "reactive_probes": ["O=CC(O)C(O)C(O)CO"],
        "alerts": "🔥 환원당 — Maillard 갈변",
    },
    "Fructose": {
        "smiles": "OCC1OC(O)(CO)C(O)C1O",
        "category": "Binder/Diluent",
        "label_ko": "과당 (Fructose)",
        "reactive_probes": ["OCC(=O)C(O)C(O)CO"],
        "alerts": "🔥 케토스 — Maillard 및 카보닐 반응",
    },
    "Trehalose": {
        "smiles": "OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O",
        "category": "Binder/Diluent",
        "label_ko": "트레할로스 (Trehalose)",
        "reactive_probes": [],
        "alerts": "비환원당 — Maillard 위험 낮음 (동결건조 안정화제)",
    },
    "Xylitol": {
        "smiles": "OCC(O)C(O)C(O)CO",
        "category": "Binder/Diluent",
        "label_ko": "자일리톨 (Xylitol)",
        "reactive_probes": [],
        "alerts": "당알코올 — 비환원, 흡습성 주의",
    },
    "Erythritol": {
        "smiles": "OCC(O)C(O)CO",
        "category": "Binder/Diluent",
        "label_ko": "에리스리톨 (Erythritol)",
        "reactive_probes": [],
        "alerts": "당알코올 — 비환원",
    },
    "Calcium_Sulfate_Dihydrate": {
        "smiles": "O.O.[Ca+2].[O-]S([O-])(=O)=O",
        "category": "Binder/Diluent",
        "label_ko": "황산칼슘 이수화물 (Calcium sulfate dihydrate)",
        "reactive_probes": [],
        "alerts": "결정수 방출 — 수분 민감 API 주의",
    },
    "Calcium_Carbonate_Precipitated": {
        "smiles": "[Ca+2].[O-]C([O-])=O",
        "category": "Binder/Diluent",
        "label_ko": "침강 탄산칼슘 (Calcium carbonate, precipitated)",
        "reactive_probes": [],
        "alerts": "알칼리 — 산 불안정 API 분해 가속",
    },
    "Magnesium_Carbonate": {
        "smiles": "[Mg+2].[O-]C([O-])=O",
        "category": "Binder/Diluent",
        "label_ko": "탄산마그네슘 (Magnesium carbonate)",
        "reactive_probes": [],
        "alerts": "알칼리 희석제 — 국소 pH 상승",
    },
    "Magnesium_Oxide": {
        "smiles": "[Mg+2].[O-2]",
        "category": "Binder/Diluent",
        "label_ko": "산화마그네슘 (Magnesium oxide)",
        "reactive_probes": [],
        "alerts": "강알칼리 — 염기성 가수분해 촉진",
    },
    "Powdered_Cellulose": {
        "smiles": "OCC1OC(O)C(O)C(O)C1O",
        "category": "Binder/Diluent",
        "label_ko": "분말 셀룰로오스 (Powdered cellulose)",
        "reactive_probes": [],
        "alerts": "수분 흡착 — 고체상 가수분해 가속",
    },
    "Silicified_MCC": {
        "smiles": "OCC1OC(O)C(O)C(O)C1O",
        "category": "Binder/Diluent",
        "label_ko": "규화 미결정셀룰로오스 (Silicified MCC)",
        "reactive_probes": [],
        "alerts": "MCC + 콜로이드 실리카 — 흡습 및 표면 흡착",
    },

    # ── Lubricant / Glidant ──────────────────────────────────────────────
    "Calcium_Stearate": {
        "smiles": "[Ca+2].CCCCCCCCCCCCCCCCCC(=O)[O-].CCCCCCCCCCCCCCCCCC(=O)[O-]",
        "category": "Lubricant/Glidant",
        "label_ko": "스테아르산칼슘 (Calcium stearate)",
        "reactive_probes": ["CCCCCCCCCCCCCCCCCC(=O)O"],
        "alerts": "유리 스테아르산 — 아민 API 아마이드화 · 용출 저하",
    },
    "Zinc_Stearate": {
        "smiles": "[Zn+2].CCCCCCCCCCCCCCCCCC(=O)[O-].CCCCCCCCCCCCCCCCCC(=O)[O-]",
        "category": "Lubricant/Glidant",
        "label_ko": "스테아르산아연 (Zinc stearate)",
        "reactive_probes": ["CCCCCCCCCCCCCCCCCC(=O)O"],
        "alerts": "Zn²⁺ — 킬레이트 형성 및 산화 촉매",
    },
    "Sodium_Stearate": {
        "smiles": "[Na+].CCCCCCCCCCCCCCCCCC(=O)[O-]",
        "category": "Lubricant/Glidant",
        "label_ko": "스테아르산나트륨 (Sodium stearate)",
        "reactive_probes": ["CCCCCCCCCCCCCCCCCC(=O)O"],
        "alerts": "알칼리성 비누 — 염기성 가수분해 촉진",
    },
    "Glyceryl_Behenate": {
        "smiles": "CCCCCCCCCCCCCCCCCCCCCC(=O)OCC(O)CO",
        "category": "Lubricant/Glidant",
        "label_ko": "글리세릴 베헤네이트 (Compritol)",
        "reactive_probes": ["CCCCCCCCCCCCCCCCCCCCCC(=O)OC", "CCCCCCCCCCCCCCCCCCCCCC(=O)O"],
        "alerts": "지방산 에스터 — 에스터 교환 및 가수분해",
    },
    "Sodium_Benzoate_Lubricant": {
        "smiles": "[Na+].[O-]C(=O)c1ccccc1",
        "category": "Lubricant/Glidant",
        "label_ko": "벤조산나트륨 (활택제 용도)",
        "reactive_probes": ["OC(=O)c1ccccc1"],
        "alerts": "산성 pH에서 벤조산 유리 — 아민 API 아마이드화",
    },

    # ── Disintegrant ─────────────────────────────────────────────────────
    "Alginic_Acid_Disintegrant": {
        "smiles": "OC1C(O)C(OC(C1O)C(=O)O)O",
        "category": "Disintegrant",
        "label_ko": "알긴산 (붕해제 용도)",
        "reactive_probes": ["OC1C(O)C(O)OC(C(=O)O)C1O"],
        "alerts": "다당류 산 — pH 저하 및 이온 교환",
    },
    "Polacrilin_Potassium": {
        "smiles": "[K+].CC(C)(C(=O)[O-])CC(C)(C)C(=O)O",
        "category": "Disintegrant",
        "label_ko": "폴라크릴린 칼륨 (Amberlite IRP-88)",
        "reactive_probes": ["CC(C)(C)C(=O)O"],
        "alerts": "이온교환 수지 — 염기성 API 이온쌍 형성 · 용출 지연",
    },

    # ── Polymer / Film coating ───────────────────────────────────────────
    "Polyethylene_Glycol_400": {
        "smiles": "OCCOCCOCCOCCO",
        "category": "Polymer/Film Coating",
        "label_ko": "폴리에틸렌글리콜 (PEG 400)",
        "reactive_probes": ["OC=O", "OO", "O=CC"],
        "alerts": "🔥 포름산·포름알데히드·과산화물 — 포르밀화 및 산화",
    },
    "Polyethylene_Glycol_6000": {
        "smiles": "OCCOCCOCCOCCO",
        "category": "Polymer/Film Coating",
        "label_ko": "폴리에틸렌글리콜 (PEG 6000)",
        "reactive_probes": ["OC=O", "OO", "O=CC"],
        "alerts": "🔥 포름산·과산화물 잔류 — NDSRI 전구체",
    },
    # ── 포비돈 등급 ──────────────────────────────────────────────────────
    #
    # ⚠️ K 값은 **분자량·점도 등급**입니다. 반응하는 자리(N-비닐피롤리돈
    #    반복 단위)는 등급과 무관하게 같으므로, 이 도구에서는 K12·K17·
    #    K25·K30·K90이 **동일한 예측**을 냅니다. 실제 처방에 쓰는 등급을
    #    그대로 고를 수 있도록 목록에만 추가한 것입니다.
    #
    #    등급에 따라 실제로 달라지는 것은 잔류 과산화물 «양»인데,
    #    이 엔진은 양을 모델링하지 않습니다.
    "Povidone_K12": {
        "smiles": "CCN1CCCC1=O",
        "category": "Polymer/Film Coating",
        "label_ko": "포비돈 K12 (PVP K12)",
        "reactive_probes": ["C(=O)O", "OO"],
        "alerts": "🔥 과산화물 잔류 — 산화 불순물 촉진 (저분자 등급)",
    },
    "Povidone_K17": {
        "smiles": "CCN1CCCC1=O",
        "category": "Polymer/Film Coating",
        "label_ko": "포비돈 K17 (PVP K17)",
        "reactive_probes": ["C(=O)O", "OO"],
        "alerts": "🔥 과산화물 잔류 — 산화 불순물 촉진",
    },
    "Povidone_K25": {
        "smiles": "CCN1CCCC1=O",
        "category": "Polymer/Film Coating",
        "label_ko": "포비돈 K25 (PVP K25)",
        "reactive_probes": ["C(=O)O", "OO"],
        "alerts": "🔥 과산화물 잔류 — 산화 불순물 촉진",
    },
    "Copovidone": {
        "smiles": "CCN1CCCC1=O.CC(=O)OC=C",
        "category": "Polymer/Film Coating",
        "label_ko": "코포비돈 (Copovidone / Kollidon VA64)",
        "reactive_probes": ["C(=O)O", "OO", "CC(=O)O"],
        "alerts": "🔥 비닐피롤리돈·비닐아세테이트 공중합체 — 과산화물 + 유리 아세트산",
    },

    "Polyvinyl_Alcohol": {
        "smiles": "OC(C)CC(C)O",
        "category": "Polymer/Film Coating",
        "label_ko": "폴리비닐알코올 (PVA)",
        "reactive_probes": ["CC(=O)O", "O=CC"],
        "alerts": "잔류 아세테이트 및 알데하이드 — 아민 API 반응",
    },
    "Methylcellulose": {
        "smiles": "COCC1OC(O)C(OC)C(O)C1O",
        "category": "Polymer/Film Coating",
        "label_ko": "메틸셀룰로오스 (Methylcellulose)",
        "reactive_probes": [],
        "alerts": "수분 보유 — 고체상 가수분해 가속",
    },
    "Hydroxyethyl_Cellulose": {
        "smiles": "OCCOCC1OC(O)C(O)C(O)C1O",
        "category": "Polymer/Film Coating",
        "label_ko": "히드록시에틸셀룰로오스 (HEC)",
        "reactive_probes": [],
        "alerts": "친수성 고분자 — 수분 활성도 증가",
    },
    "Eudragit_E": {
        "smiles": "CC(C)(C(=O)OCCN(C)C)CC(C)(C)C(=O)OC",
        "category": "Polymer/Film Coating",
        "label_ko": "아미노알킬메타크릴레이트 공중합체 E (Eudragit E)",
        "reactive_probes": ["CN(C)CCO", "CC(=O)OC"],
        "alerts": "3급 아민 함유 — 과산화물 존재 시 N-oxide, 산성 용해",
    },
    "Eudragit_RS": {
        # 2026-09-03 수정: 4급 암모늄만 있고 대이온이 없어 순전하 +1
        #   (C15H30NO4⁺)이었습니다. 암모니오메타크릴레이트 공중합체는
        #   염화물염이므로 [Cl-]를 더해 중성으로 맞춥니다.
        "smiles": "CC(C)(C(=O)OCC[N+](C)(C)C)CC(C)(C)C(=O)OC.[Cl-]",
        "category": "Polymer/Film Coating",
        "label_ko": "암모니오메타크릴레이트 공중합체 (Eudragit RS)",
        "reactive_probes": ["CC(=O)OC"],
        "alerts": "4급 암모늄 — 음이온 API 이온쌍 형성",
    },
    "Shellac": {
        "smiles": "OC(=O)CCCCCCC(O)C(O)CCCCCCC",
        "category": "Polymer/Film Coating",
        "label_ko": "셸락 (Shellac)",
        "reactive_probes": ["OC(=O)CCCCCCCC"],
        "alerts": "천연 수지산 — 유리 카복실산으로 아민 API 반응 · 경시 불용화",
    },
    "Carbomer": {
        "smiles": "OC(=O)C(C)CC(C)C(=O)O",
        "category": "Polymer/Film Coating",
        "label_ko": "카보머 (Carbopol)",
        "reactive_probes": ["CC(C)C(=O)O"],
        "alerts": "폴리아크릴산 — 염기성 API 이온쌍 · 강한 산성",
    },

    # ── Surfactant ───────────────────────────────────────────────────────
    "Polysorbate_40": {
        "smiles": "CCCCCCCCCCCCCCCC(=O)OCCOCCO",
        "category": "Surfactant",
        "label_ko": "폴리소르베이트 40 (Tween 40)",
        "reactive_probes": ["CCCCCCCCCCCCCCCC(=O)OC", "OO", "OC=O"],
        "alerts": "에스터 계면활성제 — 과산화물 및 가수분해",
    },
    "Sorbitan_Monooleate": {
        "smiles": "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(O)C1OCC(O)C1O",
        "category": "Surfactant",
        "label_ko": "소르비탄 모노올레이트 (Span 80)",
        "reactive_probes": ["CCCCCCCC/C=C\\CCCCCCCC(=O)O", "OO"],
        "alerts": "🔥 불포화 지방산 에스터 — 자동산화 과산화물",
    },
    "Docusate_Sodium": {
        "smiles": "CCCCC(CC)COC(=O)CC(C(=O)OCC(CC)CCCC)S(=O)(=O)[O-].[Na+]",
        "category": "Surfactant",
        "label_ko": "도큐세이트 나트륨 (Docusate sodium)",
        "reactive_probes": ["CCCCC(CC)COC(=O)C"],
        "alerts": "설포숙시네이트 에스터 — 가수분해 시 유리산",
    },
    "Cetrimide": {
        "smiles": "CCCCCCCCCCCCCC[N+](C)(C)C.[Br-]",
        "category": "Surfactant",
        "label_ko": "세트리마이드 (Cetrimide)",
        "reactive_probes": [],
        "alerts": "4급 암모늄 — 음이온 API 이온쌍 불순물",
    },
    "Poloxamer_188": {
        "smiles": "OCCOCC(C)OCCO",
        "category": "Surfactant",
        "label_ko": "폴록사머 188 (Pluronic F68)",
        "reactive_probes": ["OO", "OC=O"],
        "alerts": "폴리에터 — 자동산화 과산화물 및 포름산",
    },

    # ── Buffer / pH modifier ─────────────────────────────────────────────
    "Fumaric_Acid": {
        "smiles": "OC(=O)/C=C/C(=O)O",
        "category": "Buffer/pH Modifier",
        "label_ko": "푸마르산 (Fumaric acid)",
        "reactive_probes": ["OC(=O)/C=C/C(=O)O"],
        "alerts": "🔥 α,β-불포화 이산 — 아민·티올 API Michael 부가",
    },
    "Maleic_Acid": {
        "smiles": "OC(=O)/C=C\\C(=O)O",
        "category": "Buffer/pH Modifier",
        "label_ko": "말레산 (Maleic acid)",
        "reactive_probes": ["OC(=O)/C=C\\C(=O)O"],
        "alerts": "🔥 α,β-불포화 이산 — Michael 부가 · 강한 산성",
    },
    "Succinic_Acid": {
        "smiles": "OC(=O)CCC(=O)O",
        "category": "Buffer/pH Modifier",
        "label_ko": "숙신산 (Succinic acid)",
        "reactive_probes": ["OC(=O)CCC(=O)O"],
        "alerts": "이산 — 아민 API 염 형성 및 아마이드화",
    },
    "Malic_Acid": {
        "smiles": "OC(CC(=O)O)C(=O)O",
        "category": "Buffer/pH Modifier",
        "label_ko": "말산 (Malic acid)",
        "reactive_probes": ["OC(CC(=O)O)C(=O)O"],
        "alerts": "히드록시산 — 에스터화 및 아마이드화",
    },
    "Lactic_Acid": {
        "smiles": "CC(O)C(=O)O",
        "category": "Buffer/pH Modifier",
        "label_ko": "젖산 (Lactic acid)",
        "reactive_probes": ["CC(O)C(=O)O"],
        "alerts": "히드록시산 — 락타이드 형성 및 에스터화",
    },
    "Acetic_Acid_Glacial": {
        "smiles": "CC(=O)O",
        "category": "Buffer/pH Modifier",
        "label_ko": "빙초산 (Acetic acid, glacial)",
        "reactive_probes": ["CC(=O)O"],
        "alerts": "아민 API 아세틸화 · 알코올 API 에스터화",
    },
    "Sodium_Acetate": {
        "smiles": "[Na+].CC(=O)[O-]",
        "category": "Buffer/pH Modifier",
        "label_ko": "아세트산나트륨 (Sodium acetate)",
        "reactive_probes": ["CC(=O)O"],
        "alerts": "완충염 — 유리 아세트산으로 아세틸화 가능",
    },
    "Meglumine": {
        "smiles": "CNCC(O)C(O)C(O)C(O)CO",
        "category": "Buffer/pH Modifier",
        "label_ko": "메글루민 (Meglumine)",
        "reactive_probes": ["CNCC(O)C(O)C(O)C(O)CO"],
        "alerts": "🔥 2급 아민 + 폴리올 — 니트로소화 위험 및 Maillard",
    },
    "Arginine": {
        "smiles": "NC(CCCNC(N)=N)C(=O)O",
        "category": "Buffer/pH Modifier",
        "label_ko": "아르기닌 (Arginine)",
        "reactive_probes": ["NC(CCCNC(N)=N)C(=O)O"],
        "alerts": "1급 아민 + 구아니딘 — Maillard 및 아마이드화",
    },
    "Glycine": {
        "smiles": "NCC(=O)O",
        "category": "Buffer/pH Modifier",
        "label_ko": "글리신 (Glycine)",
        "reactive_probes": ["NCC(=O)O"],
        "alerts": "1급 아민 — Maillard 및 아미노산 부가",
    },

    # ── Sweetener / Flavour ──────────────────────────────────────────────
    "Acesulfame_Potassium": {
        "smiles": "[K+].CC1=CC(=O)[N-]S(=O)(=O)O1",
        "category": "Sweetener/Flavor",
        "label_ko": "아세설팜칼륨 (Acesulfame K)",
        "reactive_probes": [],
        "alerts": "설팜산 유도체 — 산성 조건에서 가수분해",
    },
    "Neotame": {
        "smiles": "CC(C)(C)CCNC(CC(=O)O)C(=O)NC(Cc1ccccc1)C(=O)OC",
        "category": "Sweetener/Flavor",
        "label_ko": "네오탐 (Neotame)",
        "reactive_probes": ["COC(=O)C"],
        "alerts": "펩타이드 에스터 — 가수분해 및 diketopiperazine",
    },
    "Sodium_Cyclamate": {
        "smiles": "[Na+].[O-]S(=O)(=O)NC1CCCCC1",
        "category": "Sweetener/Flavor",
        "label_ko": "시클라메이트 나트륨 (Sodium cyclamate)",
        "reactive_probes": [],
        "alerts": "설팜산염 — 아질산 존재 시 시클로헥실아민 유래 주의",
    },
    "Vanillin": {
        "smiles": "O=Cc1ccc(O)c(OC)c1",
        "category": "Sweetener/Flavor",
        "label_ko": "바닐린 (Vanillin)",
        "reactive_probes": ["O=Cc1ccc(O)c(OC)c1"],
        "alerts": "🔥 방향족 알데하이드 — 1급 아민 API Schiff base",
    },
    "Menthol": {
        "smiles": "CC(C)C1CCC(C)CC1O",
        "category": "Sweetener/Flavor",
        "label_ko": "멘톨 (Menthol)",
        "reactive_probes": [],
        "alerts": "테르펜 알코올 — 공융 혼합물 형성 주의",
    },

    # ── Preservative / Antioxidant ───────────────────────────────────────
    "Butylparaben": {
        "smiles": "CCCCOC(=O)c1ccc(O)cc1",
        "category": "Preservative/Antioxidant",
        "label_ko": "부틸파라벤 (Butylparaben)",
        "reactive_probes": ["COC(=O)c1ccc(O)cc1", "OC(=O)c1ccc(O)cc1"],
        "alerts": "에스터 방부제 — 아미놀리시스 및 가수분해",
    },
    "Ethylparaben": {
        "smiles": "CCOC(=O)c1ccc(O)cc1",
        "category": "Preservative/Antioxidant",
        "label_ko": "에틸파라벤 (Ethylparaben)",
        "reactive_probes": ["COC(=O)c1ccc(O)cc1", "OC(=O)c1ccc(O)cc1"],
        "alerts": "에스터 방부제 — 아미놀리시스 및 가수분해",
    },
    "Sorbic_Acid": {
        "smiles": "C/C=C/C=C/C(=O)O",
        "category": "Preservative/Antioxidant",
        "label_ko": "소르빈산 (Sorbic acid)",
        "reactive_probes": ["C/C=C/C=C/C(=O)O"],
        "alerts": "🔥 공액 불포화산 — Michael 부가 및 산화 갈변",
    },
    "Sodium_Sulfite": {
        "smiles": "[Na+].[Na+].[O-]S([O-])=O",
        "category": "Preservative/Antioxidant",
        "label_ko": "아황산나트륨 (Sodium sulfite)",
        "reactive_probes": ["OS(=O)O", "[O-]S(=O)O"],
        "alerts": "🔥 환원제 — 알데하이드·케톤 API 중아황산 부가물",
    },
    "Sodium_Bisulfite": {
        "smiles": "[Na+].OS([O-])=O",
        "category": "Preservative/Antioxidant",
        "label_ko": "아황산수소나트륨 (Sodium bisulfite)",
        "reactive_probes": ["OS(=O)O", "[O-]S(=O)O"],
        "alerts": "🔥 환원제 — 카보닐 API 부가물 · 이황화 결합 환원",
    },
    "Sodium_Thiosulfate": {
        "smiles": "[Na+].[Na+].[O-]S([O-])(=O)=S",
        "category": "Preservative/Antioxidant",
        "label_ko": "티오황산나트륨 (Sodium thiosulfate)",
        "reactive_probes": [],
        "alerts": "환원제 — 산화성 API 환원 및 이황화 절단",
    },
    "Propyl_Gallate": {
        "smiles": "CCCOC(=O)c1cc(O)c(O)c(O)c1",
        "category": "Preservative/Antioxidant",
        "label_ko": "프로필갈레이트 (Propyl gallate)",
        "reactive_probes": ["COC(=O)c1cc(O)c(O)c(O)c1", "OC(=O)c1cc(O)c(O)c(O)c1"],
        "alerts": "폴리페놀 에스터 — 금속 착화 및 산화 착색",
    },
    "Alpha_Tocopherol": {
        "smiles": "CC(C)CCCC(C)CCCC(C)CCCC1(C)CCc2c(C)c(O)c(C)c(C)c2O1",
        "category": "Preservative/Antioxidant",
        "label_ko": "α-토코페롤 (Vitamin E)",
        "reactive_probes": [],
        "alerts": "지용성 항산화제 — 산화 시 퀴논 부산물",
    },
    "Thimerosal": {
        "smiles": "CC[Hg]Sc1ccccc1C(=O)[O-].[Na+]",
        "category": "Preservative/Antioxidant",
        "label_ko": "티메로살 (Thimerosal)",
        "reactive_probes": ["OC(=O)c1ccccc1S"],
        "alerts": "유기수은 — 티올 API와 교환 반응 · 사용 제한",
    },
    "Chlorocresol": {
        "smiles": "Cc1ccc(O)c(Cl)c1",
        "category": "Preservative/Antioxidant",
        "label_ko": "클로로크레졸 (Chlorocresol)",
        "reactive_probes": [],
        "alerts": "할로겐 페놀 — 고분자 흡착 및 산화 착색",
    },

    # ── Solvent / Co-solvent ─────────────────────────────────────────────
    "Polyethylene_Glycol_300": {
        "smiles": "OCCOCCOCCO",
        "category": "Solvent/Co-solvent",
        "label_ko": "폴리에틸렌글리콜 (PEG 300)",
        "reactive_probes": ["OC=O", "OO", "O=CC"],
        "alerts": "🔥 포름산·과산화물 — 포르밀화 및 산화",
    },
    "Dimethyl_Sulfoxide": {
        "smiles": "CS(C)=O",
        "category": "Solvent/Co-solvent",
        "label_ko": "디메틸설폭사이드 (DMSO)",
        "reactive_probes": [],
        "alerts": "산성 조건 가열 시 포름알데히드 유리 가능",
    },
    "N_Methylpyrrolidone": {
        "smiles": "CN1CCCC1=O",
        "category": "Solvent/Co-solvent",
        "label_ko": "N-메틸피롤리돈 (NMP)",
        "reactive_probes": [],
        "alerts": "락탐 용매 — 가수분해 시 메틸아미노부티르산",
    },
    "Benzyl_Alcohol_Solvent": {
        "smiles": "OCc1ccccc1",
        "category": "Solvent/Co-solvent",
        "label_ko": "벤질알코올 (용매 용도)",
        "reactive_probes": ["O=Cc1ccccc1"],
        "alerts": "산화 → 벤즈알데하이드 — 1급 아민 Schiff base",
    },
    "Castor_Oil_Polyoxyl_35": {
        "smiles": "CCCCCCC(O)CC/C=C\\CCCCCCCC(=O)OCCOCCO",
        "category": "Solvent/Co-solvent",
        "label_ko": "폴리옥실 35 피마자유 (Cremophor EL)",
        "reactive_probes": ["CCCCCCC(O)CC/C=C\\CCCCCCCC(=O)O", "OO"],
        "alerts": "🔥 불포화 지방산 에스터 — 과산화물 축적",
    },
}


# ===========================================================================
# 2026-09-03 추가 — `mw_kind` 판정표
#
# 왜 필요한가 —
#   고분자·무기물의 `smiles`는 **대표 반복 단위**인데, 그 대용 SMILES도
#   RDKit에서는 멀쩡한 분자량을 반환합니다. 그래서 `p17._mole_ratio`의
#   `if not exc_mw` 가드가 한 번도 발동하지 않고 **단량체 분자량으로
#   몰비가 계산**되고 있었습니다 (PEG 4000: 보고 2.09 vs 실제 ≈0.079).
#
#   **단량체 SMILES는 반응성 탐침 용도로만 쓰며 몰비 계산에 쓰지 않습니다.**
#
# 값은 p70과 같은 세 가지뿐입니다.
#   "molecular"          단일 분자량이 정의되는 실제 화합물 — 몰비 가능
#   "monomer_surrogate"  `smiles`가 고분자의 대표 반복 단위 — 몰비 불가
#   "undefined"          광물·망상 고체·조성 미정 혼합물 — 몰비 불가
#
# 항목마다 딕셔너리에 직접 적는 대신 **한곳에 모은 판정표**를 두고
# import 시점에 각 항목에 `mw_kind` 필드로 찍습니다(아래 `_stamp_mw_kind`).
# 97건을 흩어 놓으면 누락·불일치를 눈으로 확인할 수 없기 때문입니다.
# 결과적으로 각 항목 dict는 `mw_kind` 키를 그대로 가집니다.
# ===========================================================================

MW_KIND_BY_KEY: dict[str, str] = {
    # ── STRUCTURE_PATCH 대상 (p71 FDA IID 보충 목록 + Mg 스테아레이트) ──
    # 단백질·다당류 고분자 — 반복 단위 대용
    "Gelatin": "monomer_surrogate",
    "Xanthan_Gum": "monomer_surrogate",
    "Carboxymethylcellulose_Sodium": "monomer_surrogate",
    "Alginic_Acid": "monomer_surrogate",
    "Acacia_Gum": "monomer_surrogate",
    # 무기염 — 분자량 정의됨
    "Sodium_Chloride": "molecular",
    "Calcium_Phosphate_Tribasic": "molecular",
    "Sodium_Phosphate_Dibasic": "molecular",
    "Potassium_Phosphate": "molecular",
    "Sodium_Hydroxide": "molecular",
    "Hydrochloric_Acid": "molecular",
    "Sodium_Metabisulfite": "molecular",
    "Saccharin_Sodium": "molecular",
    "Sodium_Benzoate": "molecular",
    "Potassium_Sorbate": "molecular",
    "Butylated_Hydroxyanisole_BHA": "molecular",
    "Dextrose": "molecular",
    "Isopropyl_Alcohol": "molecular",
    "Benzyl_Benzoate": "molecular",
    "Magnesium_Stearate": "molecular",
    "Glyceryl_Monostearate": "molecular",
    # TiO2 는 망상 고체
    "Titanium_Dioxide": "undefined",
    # 조성이 정해지지 않은 천연 혼합물
    "Hydrogenated_Castor_Oil": "undefined",
    "Light_Mineral_Oil": "undefined",
    "Lecithin": "undefined",
    "Peppermint_Oil": "undefined",
    # 셀룰로오스 유도체 고분자
    "Hypromellose_Phthalate": "monomer_surrogate",
    "Cellulose_Acetate_Phthalate": "monomer_surrogate",
    "Povidone_K90": "monomer_surrogate",
    "Polysorbate_60": "monomer_surrogate",

    # ── ADDITIONAL_EXCIPIENTS ────────────────────────────────────────────
    # Binder / Diluent
    "Maltodextrin": "monomer_surrogate",
    "Maltose": "molecular",
    "Fructose": "molecular",
    "Trehalose": "molecular",
    "Xylitol": "molecular",
    "Erythritol": "molecular",
    "Calcium_Sulfate_Dihydrate": "molecular",
    "Calcium_Carbonate_Precipitated": "molecular",
    "Magnesium_Carbonate": "molecular",
    "Magnesium_Oxide": "molecular",
    "Powdered_Cellulose": "monomer_surrogate",
    "Silicified_MCC": "monomer_surrogate",
    # Lubricant / Glidant
    "Calcium_Stearate": "molecular",
    "Zinc_Stearate": "molecular",
    "Sodium_Stearate": "molecular",
    "Glyceryl_Behenate": "molecular",
    "Sodium_Benzoate_Lubricant": "molecular",
    # Disintegrant
    "Alginic_Acid_Disintegrant": "monomer_surrogate",
    "Polacrilin_Potassium": "monomer_surrogate",
    # Polymer / Film coating — 전부 반복 단위 대용
    "Polyethylene_Glycol_400": "monomer_surrogate",
    "Polyethylene_Glycol_6000": "monomer_surrogate",
    "Povidone_K12": "monomer_surrogate",
    "Povidone_K17": "monomer_surrogate",
    "Povidone_K25": "monomer_surrogate",
    "Copovidone": "monomer_surrogate",
    "Polyvinyl_Alcohol": "monomer_surrogate",
    "Methylcellulose": "monomer_surrogate",
    "Hydroxyethyl_Cellulose": "monomer_surrogate",
    "Eudragit_E": "monomer_surrogate",
    "Eudragit_RS": "monomer_surrogate",
    "Carbomer": "monomer_surrogate",
    "Shellac": "undefined",           # 천연 수지산 혼합물
    # Surfactant
    "Polysorbate_40": "monomer_surrogate",
    "Sorbitan_Monooleate": "molecular",
    "Docusate_Sodium": "molecular",
    "Cetrimide": "molecular",
    "Poloxamer_188": "monomer_surrogate",
    # Buffer / pH modifier — 전부 저분자
    "Fumaric_Acid": "molecular",
    "Maleic_Acid": "molecular",
    "Succinic_Acid": "molecular",
    "Malic_Acid": "molecular",
    "Lactic_Acid": "molecular",
    "Acetic_Acid_Glacial": "molecular",
    "Sodium_Acetate": "molecular",
    "Meglumine": "molecular",
    "Arginine": "molecular",
    "Glycine": "molecular",
    # Sweetener / Flavour
    "Acesulfame_Potassium": "molecular",
    "Neotame": "molecular",
    "Sodium_Cyclamate": "molecular",
    "Vanillin": "molecular",
    "Menthol": "molecular",
    # Preservative / Antioxidant
    "Butylparaben": "molecular",
    "Ethylparaben": "molecular",
    "Sorbic_Acid": "molecular",
    "Sodium_Sulfite": "molecular",
    "Sodium_Bisulfite": "molecular",
    "Sodium_Thiosulfate": "molecular",
    "Propyl_Gallate": "molecular",
    "Alpha_Tocopherol": "molecular",
    "Thimerosal": "molecular",
    "Chlorocresol": "molecular",
    # Solvent / Co-solvent
    "Polyethylene_Glycol_300": "monomer_surrogate",
    "Dimethyl_Sulfoxide": "molecular",
    "N_Methylpyrrolidone": "molecular",
    "Benzyl_Alcohol_Solvent": "molecular",
    "Castor_Oil_Polyoxyl_35": "monomer_surrogate",
}

# 판정표에 없는 항목이 오면 «몰비 계산 금지» 쪽으로 떨어뜨립니다
MW_KIND_DEFAULT = "undefined"


def _stamp_mw_kind() -> None:
    """모듈 데이터 각 항목에 `mw_kind` 필드를 직접 찍습니다."""
    for key, patch in STRUCTURE_PATCH.items():
        patch.setdefault("mw_kind", MW_KIND_BY_KEY.get(key, MW_KIND_DEFAULT))
    for key, entry in ADDITIONAL_EXCIPIENTS.items():
        entry.setdefault("mw_kind", MW_KIND_BY_KEY.get(key, MW_KIND_DEFAULT))


_stamp_mw_kind()


# ===========================================================================
# 적용
# ===========================================================================


def apply_extensions(library: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """기존 라이브러리에 구조 보강과 신규 항목을 반영합니다.

    원본 dict를 바꾸지 않고 새 dict를 돌려줍니다.
    """
    merged: dict[str, dict[str, Any]] = {}
    for key, data in library.items():
        entry = dict(data)
        patch = STRUCTURE_PATCH.get(key)
        if patch:
            # 이미 유효한 값이 있으면 덮어쓰지 않습니다
            if not str(entry.get("smiles") or "").strip():
                entry["smiles"] = patch.get("smiles", "")
            elif not _parses(str(entry.get("smiles"))):
                entry["smiles"] = patch.get("smiles", entry.get("smiles", ""))
            if not entry.get("reactive_probes"):
                entry["reactive_probes"] = list(patch.get("reactive_probes") or [])
            entry.setdefault("source", data.get("source", "PRIME-DS Master"))
            entry["structure_source"] = "PRIME-DS extension (PubChem)"
            # 2026-09-03: 몰비 방어 메타데이터도 함께 넘깁니다.
            #   원본에 이미 값이 있으면 그것을 존중합니다(p70이 더 구체적).
            for meta_key in ("mw_kind", "repeat_unit_smiles", "structure_note"):
                if meta_key in patch and not entry.get(meta_key):
                    entry[meta_key] = patch[meta_key]

        # 프로브 보충 — 기존 것을 지우지 않고 빠진 것만 더합니다
        extra = PROBE_PATCH.get(key)
        if extra:
            probes = list(entry.get("reactive_probes") or [])
            for smi in extra:
                if smi not in probes:
                    probes.append(smi)
            entry["reactive_probes"] = probes

        merged[key] = entry

    for key, data in ADDITIONAL_EXCIPIENTS.items():
        if key in merged:
            continue
        entry = dict(data)
        entry.setdefault("source", "PRIME-DS extension")
        entry.setdefault("risk_alerts", entry.get("alerts", ""))
        merged[key] = entry

    # 2026-09-03: 어떤 경로로 들어왔든 **모든 항목**이 mw_kind를 갖게 합니다.
    #   p17이 이 필드를 읽어 몰비 계산을 방어하므로, 빠진 항목이 있으면
    #   단량체 분자량으로 몰비가 계산되던 예전 동작으로 되돌아갑니다.
    #   판정표에 없으면 안전한 쪽("undefined")으로 떨어뜨립니다.
    for key, entry in merged.items():
        if not entry.get("mw_kind"):
            entry["mw_kind"] = MW_KIND_BY_KEY.get(key, MW_KIND_DEFAULT)
        # 고분자 대용 SMILES는 반복 단위임을 필드명으로도 드러냅니다
        if entry["mw_kind"] == "monomer_surrogate" and not entry.get("repeat_unit_smiles"):
            smi = str(entry.get("smiles") or "").strip()
            if smi:
                entry["repeat_unit_smiles"] = smi

    return merged


def _parses(smiles: str) -> bool:
    try:
        from rdkit import Chem

        return Chem.MolFromSmiles(smiles) is not None
    except Exception:
        return False


def extension_stats() -> dict[str, Any]:
    return {
        "patched": len(STRUCTURE_PATCH),
        "added": len(ADDITIONAL_EXCIPIENTS),
    }


# ---------------------------------------------------------------------------
# 2026-09-03 추가 — 전수 구조 검증
#
# p19는 143종 확장 데이터의 절반 이상을 담고 있는데, 지금까지 SMILES가
# 맞는지 확인하는 코드가 없었습니다. 메타중아황산나트륨이 dithionite로,
# CMC나트륨과 Eudragit RS가 전하 불균형으로 들어와 있었습니다.
# ---------------------------------------------------------------------------


def validate_extension_smiles(*, log: bool = True) -> dict[str, Any]:
    """`STRUCTURE_PATCH` · `ADDITIONAL_EXCIPIENTS`의 구조를 전수 검사합니다.

    파싱 실패 · 순전하 불균형 · mw_kind 누락/오값을 찾습니다.
    반환 —
        {"total", "smiles_ok", "smiles_failed", "probe_failed",
         "charge_imbalanced", "bad_mw_kind"}
    """
    import logging

    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.error")
    logger = logging.getLogger(__name__)
    allowed = ("molecular", "monomer_surrogate", "undefined")

    smiles_ok = 0
    smiles_failed: list[tuple[str, str]] = []
    probe_failed: list[tuple[str, str]] = []
    charge_imbalanced: list[tuple[str, int]] = []
    bad_mw_kind: list[tuple[str, Any]] = []

    sources = (("STRUCTURE_PATCH", STRUCTURE_PATCH), ("ADDITIONAL", ADDITIONAL_EXCIPIENTS))
    total = 0
    for src_name, entries in sources:
        for key, data in entries.items():
            total += 1
            label = f"{src_name}:{key}"
            smi = str(data.get("smiles") or "").strip()
            if smi:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    smiles_failed.append((label, smi))
                    if log:
                        logger.warning("부형제 %s: smiles 파싱 실패 %r", label, smi)
                else:
                    smiles_ok += 1
                    charge = Chem.GetFormalCharge(mol)
                    if charge != 0:
                        charge_imbalanced.append((label, charge))
                        if log:
                            logger.warning(
                                "부형제 %s: 순전하 %+d — 대이온이 빠졌을 수 있습니다 (%s)",
                                label, charge, smi,
                            )
            for probe in data.get("reactive_probes") or []:
                if Chem.MolFromSmiles(str(probe)) is None:
                    probe_failed.append((label, str(probe)))
                    if log:
                        logger.warning("부형제 %s: reactive_probe 파싱 실패 %r", label, probe)

            kind = data.get("mw_kind")
            if kind not in allowed:
                bad_mw_kind.append((label, kind))
                if log:
                    logger.warning("부형제 %s: mw_kind 값이 잘못됐거나 없습니다 %r", label, kind)

    return {
        "total": total,
        "smiles_ok": smiles_ok,
        "smiles_failed": smiles_failed,
        "probe_failed": probe_failed,
        "charge_imbalanced": charge_imbalanced,
        "bad_mw_kind": bad_mw_kind,
    }
