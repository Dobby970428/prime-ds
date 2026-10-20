"""
제약 제형(DP) 업계 표준 부형제 마스터 풀 — USP/NF · Handbook of Pharmaceutical Excipients 기준.

각 항목: smiles, category, label_ko, reactive_probes, alerts, mw_kind
"""

# 부형제·완제의약품(DP)
# 파이썬으로 쓴 엑셀 표- 함수, 클래스 없음. 오직 딕셔너리
# reactive_probes: 부형제 전체 구조가 아니라 실제로 반응하는 조각을 따로 적어둔 것
# 9개 카테고리로 분류: Binder/Diluent, Lubricant/Glidant, Disintegrant, Polymer/Film Coating, Surfactant, Buffer/pH Modifier, Sweetener/Flavor, Preservative/Antioxidant, Solvent/Co-solvent.

# ===========================================================================
# ⚠️ 2026-09-03 추가 — `mw_kind`: 분자량을 몰비 계산에 써도 되는지 표시
#
# **단량체 SMILES는 반응성 탐침(reactivity probe) 용도로만 쓰며,
#   몰비(mole ratio) 계산에 쓰지 않습니다.**
#
# 고분자·무기물은 전체 구조를 SMILES로 쓸 수 없어 **대표 반복 단위**를
# 적어 두었습니다. 그런데 그 대용 SMILES도 RDKit에서는 멀쩡한 분자량을
# 반환하기 때문에, `p17._mole_ratio`의 `if not exc_mw` 가드가 한 번도
# 발동하지 않고 **단량체 분자량으로 몰비가 계산되고 있었습니다.**
#
#     PEG 4000  보고 몰비 2.09  ← 트리에틸렌글리콜(MW 150.17) 기준
#               실제 고분자 기준으로는 ≈0.079 (26배 오차)
#     크로스포비돈·MCC          애초에 몰수가 정의되지 않음
#
# 그래서 항목마다 다음 셋 중 하나를 명시합니다.
#
#     "molecular"          단일 분자량이 정의되는 실제 화합물.
#                          몰비 계산에 **써도 됩니다.**
#     "monomer_surrogate"  `smiles`가 고분자의 대표 반복 단위/조각입니다.
#                          반응 탐침으로만 쓰고 몰비에 **쓰면 안 됩니다.**
#                          실제 반복 단위는 `repeat_unit_smiles`에 둡니다.
#     "undefined"          광물·망상 고체·조성이 정해지지 않은 혼합물이라
#                          분자량 자체가 정의되지 않습니다. 몰비 **불가**.
#
# p17이 이 필드를 읽어 방어합니다. 값은 반드시 위 세 문자열 중 하나여야
# 합니다. 기존 `smiles` 키는 하위 호환을 위해 그대로 둡니다.
# ===========================================================================

from __future__ import annotations

from typing import Any

# mw_kind 허용값 — p17·p19가 함께 씁니다
MW_KIND_MOLECULAR = "molecular"
MW_KIND_MONOMER_SURROGATE = "monomer_surrogate"
MW_KIND_UNDEFINED = "undefined"
MW_KIND_VALUES: tuple[str, ...] = (
    MW_KIND_MOLECULAR,
    MW_KIND_MONOMER_SURROGATE,
    MW_KIND_UNDEFINED,
)

# 카테고리 표시 순서 (UI)
EXCIPIENT_CATEGORY_ORDER: tuple[str, ...] = (
    "Binder/Diluent",
    "Lubricant/Glidant",
    "Disintegrant",
    "Polymer/Film Coating",
    "Surfactant",
    "Buffer/pH Modifier",
    "Sweetener/Flavor",
    "Preservative/Antioxidant",
    "Solvent/Co-solvent",
)

COMPLETE_EXCIPIENT_MASTER: dict[str, dict[str, Any]] = {
    # ── 1. Binders / Diluents ──────────────────────────────────────────────
    # 2026-09-03 수정: 일수화물과 무수물이 **완전히 같은 SMILES**여서
    #   구분되지 않았습니다. 일수화물에 결정수 `.O`를 붙여 실제 조성
    #   (C12H24O12, MW 360.31)을 반영했습니다. 무수물은 C12H22O11 (342.30).
    "Lactose_Monohydrate": {
        "smiles": "OCC1OC(OC2C(O)C(O)OC(CO)C2O)C(O)C(O)C1O.O",
        "category": "Binder/Diluent",
        "label_ko": "유당 일수화물 (Lactose monohydrate)",
        "reactive_probes": ["O=CC(O)CO"],
        "alerts": "🔥 환원당 — 1·2급 아민 API Maillard 갈변/축합 불순물",
        "mw_kind": "molecular",
    },
    "Lactose_Anhydrous": {
        "smiles": "OCC1OC(OC2C(O)C(O)OC(CO)C2O)C(O)C(O)C1O",
        "category": "Binder/Diluent",
        "label_ko": "무수유당 (Lactose anhydrous)",
        "reactive_probes": ["O=CC(O)CO"],
        "alerts": "🔥 저수분 제형에서도 환원당 개방형 알데하이드 잔류 Maillard 리스크",
        "mw_kind": "molecular",
    },
    "Microcrystalline_Cellulose": {
        # 2026-09-03: 셀룰로오스 고분자 — smiles는 β-D-글루코스 반복 단위 대용입니다
        "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O",
        "repeat_unit_smiles": "C(C1C(C(C(C(O1)O)O)O)O)O",
        "category": "Binder/Diluent",
        "label_ko": "미결정 셀룰로오스 (MCC / Avicel)",
        "reactive_probes": [],
        "alerts": "수분 흡착 → 고체상 가수분해(특히 에스터·아미드 API) 가속",
        "mw_kind": "monomer_surrogate",
        "structure_note": "β-1,4-글루칸 고분자. 반복 단위 MW(180)는 몰비 계산에 쓸 수 없습니다.",
    },
    # 2026-09-03 수정: 만니톨과 소르비톨이 **입체 없이 같은 SMILES**로
    #   코딩되어 구분되지 않았습니다. 둘은 C4 부분입체이성질체입니다.
    #   D-만니톨 InChIKey FBPFZTCFMRRESA-KVTDHHQDSA-N
    #   D-소르비톨 InChIKey FBPFZTCFMRRESA-JGWLITMVSA-N (실측 확인)
    "Mannitol": {
        "smiles": "OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
        "category": "Binder/Diluent",
        "label_ko": "만니톨 (D-Mannitol)",
        "reactive_probes": ["O=CC(O)C(O)C(O)CO"],
        "alerts": "알코올 탈수 부산물 · 교차 혼입 시 산화 촉매 잔류 주의",
        "mw_kind": "molecular",
    },
    "Sorbitol": {
        "smiles": "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
        "category": "Binder/Diluent",
        "label_ko": "소르비톨 (D-Sorbitol / D-glucitol)",
        "reactive_probes": ["O=CC(O)C(O)C(O)CO"],
        "alerts": "히그로스코픽 당알코올 — Maillard 및 미량 알데하이드 형성",
        "mw_kind": "molecular",
    },
    # 2026-09-03 수정: Ca²⁺ 2개에 인산수소 1개라 전하가 +2 남아 있었고
    #   분자식이 HCa2O4P⁺²로 나왔습니다. DCP는 CaHPO4 (중성, MW 136.06).
    "Dibasic_Calcium_Phosphate": {
        "smiles": "OP(=O)([O-])[O-].[Ca+2]",
        "category": "Binder/Diluent",
        "label_ko": "인산수소칼슘 (DCP / Emcompress)",
        "reactive_probes": [],
        "alerts": "알칼리성 희석제 — 국소 pH 상승에 따른 염기성 가수분해 촉진",
        "mw_kind": "molecular",
    },
    "Calcium_Carbonate": {
        "smiles": "O=C([O-])[O-].[Ca+2]",
        "category": "Binder/Diluent",
        "label_ko": "탄산칼슘 (Calcium carbonate)",
        "reactive_probes": [],
        "alerts": "제형 내 CO₂/알칼리 유리 — 산 불안정 API 분해 가속",
        "mw_kind": "molecular",
    },
    # 2026-09-03: 옥수수 전분과 호화 전분은 **같은 α-글루칸 고분자**이고
    #   차이는 물리적 호화 공정입니다. 구조로는 구분되지 않으므로 같은
    #   반복 단위를 쓰되 `mw_kind`/`structure_note`로 명시합니다.
    "Corn_Starch": {
        "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "repeat_unit_smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "category": "Binder/Diluent",
        "label_ko": "옥수수 전분 (Corn starch)",
        "reactive_probes": ["O=CC(O)CO"],
        "alerts": "환원당 잔류 — 아민 API Maillard 및 색변 리스크",
        "mw_kind": "monomer_surrogate",
        "structure_note": "α-1,4/1,6-글루칸(아밀로스·아밀로펙틴). 반복 단위는 α-D-글루코스.",
    },
    "Pregelatinized_Starch": {
        "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "repeat_unit_smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "category": "Binder/Diluent",
        "label_ko": "호화 전분 (Pregelatinized starch)",
        "reactive_probes": ["O=CC(O)CO"],
        "alerts": "고수분 결합제 — Maillard 및 고체상 가수분해 촉매",
        "mw_kind": "monomer_surrogate",
        "structure_note": (
            "옥수수 전분과 **화학 구조가 같습니다** — 차이는 호화(물리적 전처리)이며 "
            "구조로는 구분되지 않습니다. 환원 말단이 더 많이 노출됩니다."
        ),
    },
    "Sucrose": {
        "smiles": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2OC(CO)C(O)C(O)C2O",
        "category": "Binder/Diluent",
        "label_ko": "자당 (Sucrose)",
        "reactive_probes": ["O=CC(O)CO", "O=CO"],
        "alerts": "🔥 가역적 환원당 — 열/산 분해 후 Maillard · caramelization",
        "mw_kind": "molecular",
    },
    "Isomalt": {
        # ⚠️ 2026-09-03 확인: 이 SMILES는 헵티톨(C7H16O7, MW 212.20)이며
        #    실제 이소말트(GPS/GPM 등몰 혼합물, 각 C12H24O11 MW 344.31)가
        #    아닙니다. 조성이 정해진 «혼합물»이라 단일 구조로 대표시키기
        #    어려워 구조는 그대로 두고 mw_kind를 undefined로 표시합니다.
        "smiles": "OCC(O)C(O)C(O)C(O)C(O)CO",
        "category": "Binder/Diluent",
        "label_ko": "이소말트 (Isomalt)",
        "reactive_probes": ["O=CC(O)CO"],
        "alerts": "당알코올 혼합물 — 저수분에서도 환원성 잔류물 주의",
        "mw_kind": "undefined",
        "structure_note": "GPS/GPM 등몰 혼합물. 등재 SMILES는 폴리올 대용 조각입니다.",
    },
    # ── 2. Lubricants / Glidants ───────────────────────────────────────────
    # 2026-09-03 수정: 두 번째 스테아레이트 사슬 **전체가 대괄호로 감싸져**
    #   있어 원자 스펙으로 해석됐습니다 (`[CCCCCCCCCCCCCCCCCC(=O)[O-]]`).
    #   마스터 47종 중 유일하게 MolFromSmiles가 None을 냈고,
    #   p71.probe_smiles_for_entry가 None을 조용히 버려 평가에서 사라졌습니다.
    #   바깥 대괄호 제거 → C36H70MgO4, MW 591.26 (실측 확인).
    "Magnesium_Stearate": {
        "smiles": "CCCCCCCCCCCCCCCCCC(=O)[O-].CCCCCCCCCCCCCCCCCC(=O)[O-].[Mg+2]",
        "category": "Lubricant/Glidant",
        "label_ko": "마그네슘 스테아레이트",
        "reactive_probes": ["CCCCCCCCCCCCCCCCCC(=O)O"],
        "alerts": "🔥 유리 지방산 + Mg²⁺ — 산화·가수분해 촉매, 용출 저하",
        "mw_kind": "molecular",
    },
    "Stearic_Acid": {
        "smiles": "CCCCCCCCCCCCCCCCCC(=O)O",
        "category": "Lubricant/Glidant",
        "label_ko": "스테아르산 (Stearic acid)",
        "reactive_probes": ["CCCCCCCCCCCCCCCCCC(=O)O"],
        "alerts": "아민/알코올 API와 에스테르화·염 형성 2차 불순물",
        "mw_kind": "molecular",
    },
    "Talc": {
        # ⚠️ 실제 활석은 Mg3Si4O10(OH)2 층상 광물입니다. 등재 SMILES는
        #    Mg–규산염을 흉내 낸 조각(MgO5Si2)이라 분자량이 의미 없습니다.
        "smiles": "[Mg+2].[O-][Si](=O)O[Si](=O)[O-]",
        "category": "Lubricant/Glidant",
        "label_ko": "활석 (Talc)",
        "reactive_probes": [],
        "alerts": "전이금속 trace — 라디칼 자동산화 개시 위험",
        "mw_kind": "undefined",
        "structure_note": "층상 규산염 광물 Mg3Si4O10(OH)2 — 분자량 정의 없음.",
    },
    "Colloidal_Silicon_Dioxide": {
        "smiles": "O=[Si]=O",
        "category": "Lubricant/Glidant",
        "label_ko": "콜로이드성 이산화규소 (Aerosil)",
        "reactive_probes": [],
        "alerts": "고표면적 흡착 — 수분·산화제 농축에 따른 API 불안정",
        "mw_kind": "undefined",
        "structure_note": "SiO₂ 망상 고체 — 분자량 정의 없음.",
    },
    # 2026-09-03 수정: 기존 SMILES는 스테아**로일** 엔올 에스터
    #   (`C(=O)OC=C…`, C21H37NaO4)로, 경고문이 말하는 «푸마르산 에스터»가
    #   아니었습니다. SSF는 푸마르산의 **스테아릴(옥타데실) 모노에스터**
    #   나트륨염입니다 → C22H39NaO4, MW 390.54 (실측 확인).
    "Sodium_Stearyl_Fumarate": {
        "smiles": "CCCCCCCCCCCCCCCCCCOC(=O)/C=C/C(=O)[O-].[Na+]",
        "category": "Lubricant/Glidant",
        "label_ko": "스테아릴푸마르산 나트륨 (SSF / Pruv)",
        "reactive_probes": ["OC(=O)/C=C/C(=O)O"],
        "alerts": "푸마르산 에스터 — 아민 API와 Michael-type 가역 반응 가능",
        "mw_kind": "molecular",
    },
    # ── 3. Disintegrants ───────────────────────────────────────────────────
    # 2026-09-03: 크로스포비돈은 가교 PVP — 반복 단위는 포비돈 K30과 같습니다.
    #   구조로는 구분되지 않고 «가교되어 불용»이라는 점이 다릅니다.
    "Crospovidone": {
        "smiles": "CCN1CCCC1=O",
        "repeat_unit_smiles": "CCN1CCCC1=O",
        "category": "Disintegrant",
        "label_ko": "크로스포비돈 (Crospovidone / Kollidon CL)",
        "reactive_probes": ["C(=O)O", "OO"],
        "alerts": "🔥 잔류 과산화물 — API N-oxide / sulfoxide 산화",
        "mw_kind": "monomer_surrogate",
        "structure_note": (
            "가교 폴리(N-비닐피롤리돈) — 불용성이라 분자량이 정의되지 않습니다. "
            "반복 단위는 포비돈 K30과 동일하며 차이는 가교 여부입니다."
        ),
    },
    # 2026-09-03 수정: 크로스카멜로오스와 SSG가 **둘 다 젖산나트륨**
    #   (`CC(C(=O)[O-])O.[Na+]`, C3H5NaO3)으로 코딩되어 있었습니다.
    #   셀룰로오스·전분 골격도 카복시메틸기도 없어, 구분은커녕 화학적으로
    #   틀린 구조였습니다. 카복시메틸화 글루코스 반복 단위로 교체합니다.
    #     크로스카멜로오스 → β-글루칸(셀룰로오스) 골격
    #     SSG             → α-글루칸(전분) 골격
    #   둘 다 C8H13NaO8 (MW 260.17)이지만 아노머 배치가 달라 구분됩니다.
    "Croscarmellose_Sodium": {
        "smiles": "[Na+].OC[C@H]1O[C@@H](O)[C@H](OCC(=O)[O-])[C@@H](O)[C@@H]1O",
        "repeat_unit_smiles": "[Na+].OC[C@H]1O[C@@H](O)[C@H](OCC(=O)[O-])[C@@H](O)[C@@H]1O",
        "category": "Disintegrant",
        "label_ko": "크로스카멜로오스 나트륨 (Ac-Di-Sol)",
        "reactive_probes": ["OCC(=O)O"],
        "alerts": "급속 팽윤 → 국소 pH 변동 · ambient hydrolysis 가속",
        "mw_kind": "monomer_surrogate",
        "structure_note": "가교 카복시메틸셀룰로오스 나트륨 — β-1,4-글루칸 골격의 CM 유도체.",
    },
    "Sodium_Starch_Glycolate": {
        "smiles": "[Na+].OC[C@H]1O[C@H](O)[C@H](OCC(=O)[O-])[C@@H](O)[C@@H]1O",
        "repeat_unit_smiles": "[Na+].OC[C@H]1O[C@H](O)[C@H](OCC(=O)[O-])[C@@H](O)[C@@H]1O",
        "category": "Disintegrant",
        "label_ko": "나트륨 녹말글리콜산 (SSG / Explotab)",
        "reactive_probes": ["O=CC(O)CO", "OCC(=O)O"],
        "alerts": "전분 유래 환원당 + 카르복실산 — Maillard 교차 반응",
        "mw_kind": "monomer_surrogate",
        "structure_note": "카복시메틸 전분 나트륨 — α-글루칸 골격의 CM 유도체 (셀룰로오스가 아님).",
    },
    "Low_Substituted_HPC": {
        # ⚠️ 등재 SMILES는 이소프로판올 조각으로, 셀룰로오스 골격이 없는
        #    «히드록시프로필기» 대용입니다. 반응 탐침으로만 씁니다.
        "smiles": "CC(C)O",
        "repeat_unit_smiles": "CC(C)O",
        "category": "Disintegrant",
        "label_ko": "저치환 히드록시프로필셀룰로오스 (L-HPC)",
        "reactive_probes": [],
        "alerts": "수분 흡수 촉진 — 고체상 가수분해 가속",
        "mw_kind": "monomer_surrogate",
        "structure_note": "셀룰로오스 에터 고분자. 등재 조각은 히드록시프로필 치환기 대용입니다.",
    },
    # ── 4. Polymers / Film coating ─────────────────────────────────────────
    # 2026-09-03: 아래 고분자 전부 `smiles`가 대표 반복 단위/조각입니다.
    #   PEG 4000의 «몰비 2.09»는 트리에틸렌글리콜(MW 150.17) 기준이었고
    #   실제 고분자 기준으로는 ≈0.079 — 26배 오차였습니다. mw_kind로 차단합니다.
    "HPMC": {
        "smiles": "COCCOC",
        "repeat_unit_smiles": "COCCOC",
        "category": "Polymer/Film Coating",
        "label_ko": "히드록시프로필메틸셀룰로오스 (HPMC)",
        "reactive_probes": ["C(=O)O"],
        "alerts": "막코팅 중 수분 장벽 — 코어-막 pH 불균형 분해",
        "mw_kind": "monomer_surrogate",
        "structure_note": "셀룰로오스 에터. 등재 조각은 메톡시/히드록시프로폭시 치환기 대용입니다.",
    },
    "Povidone_K30": {
        "smiles": "CCN1CCCC1=O",
        "repeat_unit_smiles": "CCN1CCCC1=O",
        "category": "Polymer/Film Coating",
        "label_ko": "포비돈 K30 (PVP K30)",
        "reactive_probes": ["C(=O)O", "OO"],
        "alerts": "🔥 과산화물 잔류 — 산화 불순물 촉진",
        "mw_kind": "monomer_surrogate",
        "structure_note": (
            "선형 폴리(N-비닐피롤리돈). K값은 분자량·점도 등급이며 반복 단위는 "
            "크로스포비돈과 같습니다(가교 여부가 다름)."
        ),
    },
    "Polyethylene_Glycol_4000": {
        "smiles": "OCCOCCOCCO",
        "repeat_unit_smiles": "OCCOCCOCCO",
        "category": "Polymer/Film Coating",
        "label_ko": "폴리에틸렌글리콜 (PEG 4000)",
        "reactive_probes": ["C(=O)O", "O=CO"],
        "alerts": "🔥 포름산/포름알데히드 — Formylation · NDSRI 전구체",
        "mw_kind": "monomer_surrogate",
        "structure_note": (
            "공칭 평균분자량 ≈4000. 등재 조각은 트리에틸렌글리콜(MW 150.17)이므로 "
            "몰비에 쓰면 약 26배 과대평가됩니다."
        ),
    },
    "Ethyl_Cellulose": {
        "smiles": "CCOCCOC",
        "repeat_unit_smiles": "CCOCCOC",
        "category": "Polymer/Film Coating",
        "label_ko": "에틸셀룰로오스 (Ethyl cellulose)",
        "reactive_probes": ["C(=O)O"],
        "alerts": "유기용매 잔류 + 산화 분해 유기산 — API 상호작용",
        "mw_kind": "monomer_surrogate",
        "structure_note": "셀룰로오스 에터. 등재 조각은 에톡시 치환기 대용입니다.",
    },
    "Hydroxypropyl_Cellulose": {
        "smiles": "CC(C)OCCOC",
        "repeat_unit_smiles": "CC(C)OCCOC",
        "category": "Polymer/Film Coating",
        "label_ko": "히드록시프로필셀룰로오스 (HPC / Klucel)",
        "reactive_probes": [],
        "alerts": "결합제·코팅 — 수분 흡착에 따른 고체상 불안정",
        "mw_kind": "monomer_surrogate",
        "structure_note": "셀룰로오스 에터. 등재 조각은 히드록시프로필 치환기 대용입니다.",
    },
    "Methacrylic_Acid_Copolymer": {
        "smiles": "CC(=C)C(=O)O",
        "repeat_unit_smiles": "CC(=C)C(=O)O",
        "category": "Polymer/Film Coating",
        "label_ko": "메타아크릴산 공중합체 (Eudragit L/S)",
        "reactive_probes": ["CC(=O)O"],
        "alerts": "pH 의존 용출 — 장내 알칼리에서 API 염기성 가수분해",
        "mw_kind": "monomer_surrogate",
        "structure_note": "메타크릴산/메타크릴산메틸 공중합체. 등재 구조는 단량체입니다.",
    },
    "Poloxamer_407": {
        "smiles": "COCCOCCOC",
        "repeat_unit_smiles": "COCCOCCOC",
        "category": "Polymer/Film Coating",
        "label_ko": "폴록사머 407 (Pluronic F127)",
        "reactive_probes": ["C(=O)O", "OO"],
        "alerts": "온도 민감 겔 — 자동산화 과산화물 잔류",
        "mw_kind": "monomer_surrogate",
        "structure_note": "PEO–PPO–PEO 블록 공중합체(공칭 ≈12600). 등재 조각은 폴리에터 대용입니다.",
    },
    # ── 5. Surfactants ─────────────────────────────────────────────────────
    "Polysorbate_80": {
        # 등재 SMILES는 폴리소르베이트가 아니라 «불포화 지방산 사슬» 조각입니다
        "smiles": "CCCCC=CCCCCCCCC(=O)O",
        "repeat_unit_smiles": "CCCCC=CCCCCCCCC(=O)O",
        "category": "Surfactant",
        "label_ko": "폴리소르베이트 80 (Tween 80)",
        "reactive_probes": ["C(=O)O", "OO"],
        "alerts": "🔥 불포화 지방산 — 과산화물 · 에폭시화 분해",
        "mw_kind": "monomer_surrogate",
        "structure_note": "폴리옥시에틸렌 소르비탄 모노올레이트(공칭 ≈1310). 등재 조각은 지방산 대용입니다.",
    },
    "Polysorbate_20": {
        "smiles": "CCCCCCCCCC(=O)O",
        "repeat_unit_smiles": "CCCCCCCCCC(=O)O",
        "category": "Surfactant",
        "label_ko": "폴리소르베이트 20 (Tween 20)",
        "reactive_probes": ["C(=O)O", "OO"],
        "alerts": "에스터 계면활성제 — 가수분해·과산화물 생성",
        "mw_kind": "monomer_surrogate",
        "structure_note": "폴리옥시에틸렌 소르비탄 모노라우레이트(공칭 ≈1228). 등재 조각은 지방산 대용입니다.",
    },
    "Sodium_Lauryl_Sulfate": {
        "smiles": "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]",
        "category": "Surfactant",
        "label_ko": "라우릴 황산나트륨 (SLS)",
        "reactive_probes": [],
        "alerts": "강한 계면활성 — 수분 접근성 증가 → 가수분해 가속",
        "mw_kind": "molecular",
    },
    "Benzalkonium_Chloride": {
        # ⚠️ 2026-09-03 확인: 등재 SMILES에 **벤질기가 없습니다**
        #    (도데실트리메틸암모늄, C15H34ClN). 실제 BAC는 알킬디메틸
        #    **벤질**암모늄 클로라이드의 동족체 혼합물(C12~C16)입니다.
        #    동족체 혼합물이라 단일 분자량이 정의되지 않습니다.
        "smiles": "CCCCCCCCCCCC[N+](C)(C)C.[Cl-]",
        "category": "Surfactant",
        "label_ko": "벤잘코늄염화물 (BAC)",
        "reactive_probes": [],
        "alerts": "4급 암모늄 — 음이온 API와 이온쌍 불순물",
        "mw_kind": "undefined",
        "structure_note": "알킬(C12~C16)디메틸벤질암모늄 클로라이드 동족체 혼합물 — 단일 분자량 없음.",
    },
    # ── 6. Buffers / pH modifiers ──────────────────────────────────────────
    "Citric_Acid": {
        "smiles": "OC(=O)CC(O)(CC(=O)O)C(=O)O",
        "category": "Buffer/pH Modifier",
        "label_ko": "구연산 (Citric acid)",
        "reactive_probes": ["OC(=O)CC(O)(CC(=O)O)C(=O)O"],
        "alerts": "산성 pH — 산 가수분해(에스터·아미드·β-락탐) 촉진",
        "mw_kind": "molecular",
    },
    "Sodium_Citrate": {
        "smiles": "C(C(=O)[O-])C(CC(=O)[O-])(C(=O)[O-])O.[Na+].[Na+].[Na+]",
        "category": "Buffer/pH Modifier",
        "label_ko": "구연산나트륨 (Sodium citrate)",
        "reactive_probes": ["OC(=O)CC(O)(CC(=O)O)C(=O)O"],
        "alerts": "완충 알칼리 — pH 6~7에서 염기성 가수분해 경로 활성",
        "mw_kind": "molecular",
    },
    "Tartaric_Acid": {
        "smiles": "OC(C(=O)O)C(O)C(=O)O",
        "category": "Buffer/pH Modifier",
        "label_ko": "타르타르산 (Tartaric acid)",
        "reactive_probes": ["OC(C(=O)O)C(O)C(=O)O"],
        "alerts": "이황기산 — 아민 API 염 형성 · 고체상 반응",
        "mw_kind": "molecular",
    },
    # 2026-09-03 수정: 탄산 이온(CO3²⁻) + Na⁺ 2개, 즉 **탄산나트륨**
    #   (CNa2O3, MW 105.99)이 코딩되어 있었습니다. 탄산수소나트륨은
    #   HCO3⁻ + Na⁺ 하나 → CHNaO3, MW 84.01 (실측 확인).
    "Sodium_Bicarbonate": {
        "smiles": "OC(=O)[O-].[Na+]",
        "category": "Buffer/pH Modifier",
        "label_ko": "탄산수소나트륨 (Sodium bicarbonate)",
        "reactive_probes": ["C(=O)O"],
        "alerts": "CO₂ 방출 — 산 불안정 API · pH 국소 상승",
        "mw_kind": "molecular",
    },
    "Tromethamine_TRIS": {
        "smiles": "NC(CO)(CO)CO",
        "category": "Buffer/pH Modifier",
        "label_ko": "트로메타민 (Tris / THAM)",
        "reactive_probes": ["NC(CO)(CO)CO"],
        "alerts": "1·2급 아민 완충제 — 알데하이드 API와 Schiff base 형성",
        "mw_kind": "molecular",
    },
    # ── 7. Sweeteners / Flavors ────────────────────────────────────────────
    "Aspartame": {
        "smiles": "COC(=O)C(N)Cc1ccccc1",
        "category": "Sweetener/Flavor",
        "label_ko": "아스파탐 (Aspartame)",
        "reactive_probes": ["COC(=O)C(N)Cc1ccccc1"],
        "alerts": "🔥 펩타이드 결합 — pH·수분에서 Asp-Phe diketopiperazine",
        "mw_kind": "molecular",
    },
    # 2026-09-03 수정: 사이클로헥산 고리에 Cl 4개짜리 구조(C6H8Cl4O2)로
    #   코딩되어 있었습니다 — 이당류도 아니고 Cl 개수도 틀렸습니다.
    #   수크랄로스는 Cl 3개의 삼염화 이당류 → C12H19Cl3O8, MW 397.64
    #   (InChIKey BAQAVOSOZGMPRM-UGDNZRGBSA-N, 실측 확인).
    "Sucralose": {
        "smiles": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O[C@]2([C@H]([C@@H]([C@H](O2)CCl)O)O)CCl)O)O)Cl)O",
        "category": "Sweetener/Flavor",
        "label_ko": "수크랄로스 (Sucralose)",
        "reactive_probes": [],
        "alerts": "염소 치환당 — 열 분해 시 HCl 유리 → 산 가수분해",
        "mw_kind": "molecular",
    },
    # ── 8. Preservatives / Antioxidants ────────────────────────────────────
    "Methylparaben": {
        "smiles": "COC(=O)c1ccc(O)cc1",
        "category": "Preservative/Antioxidant",
        "label_ko": "메틸파라벤 (Methylparaben)",
        "reactive_probes": ["COC(=O)c1ccc(O)cc1"],
        "alerts": "에스터 방부제 — 가수분해 후 p-히드록시벤조산 잔류",
        "mw_kind": "molecular",
    },
    "Propylparaben": {
        "smiles": "CCCOC(=O)c1ccc(O)cc1",
        "category": "Preservative/Antioxidant",
        "label_ko": "프로필파라벤 (Propylparaben)",
        "reactive_probes": ["CCCOC(=O)c1ccc(O)cc1"],
        "alerts": "페놀 에스터 — 알칼리 가수분해 · 페놀 유리",
        "mw_kind": "molecular",
    },
    "Benzyl_Alcohol": {
        "smiles": "OCc1ccccc1",
        "category": "Preservative/Antioxidant",
        "label_ko": "벤질알코올 (Benzyl alcohol)",
        "reactive_probes": ["O=Cc1ccccc1"],
        "alerts": "산화 → 벤즈알데하이드 — Maillard · Schiff base",
        "mw_kind": "molecular",
    },
    "Butylated_Hydroxytoluene_BHT": {
        "smiles": "Cc1cc(C(C)(C)C)c(O)c(C)c1",
        "category": "Preservative/Antioxidant",
        "label_ko": "BHT (Butylated hydroxytoluene)",
        "reactive_probes": [],
        "alerts": "항산화제 — 산화 억제하나 quinone류 분해 부산물 추적",
        "mw_kind": "molecular",
    },
    # 2026-09-03 수정: 등재 구조는 **락톤 고리도 엔다이올도 없는** 열린
    #   폴리하이드록시 케톤(C5H10O5, MW 150.13)이었습니다. 아스코르브산의
    #   환원성은 γ-락톤 위의 엔다이올(C2-C3 =C(O)-C(O)) 에서 나오므로,
    #   이 구조로는 산화·Maillard 보조 경로가 생성되지 않습니다.
    #   L-아스코르브산 C6H8O6, MW 176.12
    #   (InChIKey CIWBSHSKHKDKBQ-JLAZNSOCSA-N, 실측 확인).
    "Ascorbic_Acid": {
        "smiles": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
        "category": "Preservative/Antioxidant",
        "label_ko": "비타민 C (L-Ascorbic acid)",
        "reactive_probes": ["OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", "O=CC(O)C(=O)CO"],
        "alerts": "환원제 — 금속 촉매 하 API 산화환원 · Maillard 보조",
        "mw_kind": "molecular",
    },
    "Edetate_Disodium_EDTA": {
        # 등재 구조는 유리산(EDTA, C10H16N2O8)입니다. 이나트륨염이 아닙니다.
        "smiles": "O=C(O)CN(CCN(CC(=O)O)CC(=O)O)CC(=O)O",
        "category": "Preservative/Antioxidant",
        "label_ko": "EDTA-2Na (Edetate disodium)",
        "reactive_probes": [],
        "alerts": "금속 킬레이트 — 촉매 금속 제거/이동에 따른 분해 경로 변화",
        "mw_kind": "molecular",
        "structure_note": "등재 구조는 EDTA 유리산 — 실제 품목은 이나트륨염(MW 336.21, 이수화물 372.24).",
    },
    # ── 9. Solvents / Co-solvents ──────────────────────────────────────────
    "Propylene_Glycol": {
        "smiles": "CC(O)CO",
        "category": "Solvent/Co-solvent",
        "label_ko": "프로필렌글리콜 (PG)",
        "reactive_probes": ["CC(=O)CO"],
        "alerts": "산화 → 젖산·아세톤 — 국소 pH 변동",
        "mw_kind": "molecular",
    },
    "Glycerol": {
        "smiles": "OCC(O)CO",
        "category": "Solvent/Co-solvent",
        "label_ko": "글리세롤 (Glycerin)",
        "reactive_probes": ["O=CC(O)CO"],
        "alerts": "히그로스코픽 — 수분 활성도 증가 · Maillard 보조",
        "mw_kind": "molecular",
    },
    "Ethanol": {
        "smiles": "CCO",
        "category": "Solvent/Co-solvent",
        "label_ko": "에탄올 (Ethanol)",
        "reactive_probes": ["CC=O"],
        "alerts": "산화 → 아세트알데하이드 — 아민 API 축합",
        "mw_kind": "molecular",
    },
}

# ---------------------------------------------------------------------------
# 2026-09-03 추가 — 로드 시점 구조 검증
#
# 마그네슘 스테아레이트의 SMILES가 문법 오류였는데도 **전 계층이 무음**
# 이었습니다. `p71.probe_smiles_for_entry`는 MolFromSmiles가 None이면 그냥
# `return` 하고, 깨진 문자열은 카탈로그·부형제 보고서에 그대로 출력됐습니다.
# 이제 라이브러리를 쓰기 전에 전 항목을 파싱해 실패를 남깁니다.
# ---------------------------------------------------------------------------

def validate_master_library(
    library: dict[str, dict[str, Any]] | None = None,
    *,
    log: bool = True,
) -> dict[str, Any]:
    """전 항목의 `smiles`·`reactive_probes`·`mw_kind`를 검증합니다.

    반환 —
        {"total", "smiles_ok", "smiles_failed", "probe_failed",
         "bad_mw_kind", "missing_mw_kind"}
    실패 항목은 `logging.warning`으로도 남깁니다 (`log=False`면 생략).
    """
    import logging

    from rdkit import Chem
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.error")

    lib = COMPLETE_EXCIPIENT_MASTER if library is None else library
    logger = logging.getLogger(__name__)

    smiles_ok = 0
    smiles_failed: list[tuple[str, str]] = []
    probe_failed: list[tuple[str, str]] = []
    bad_mw_kind: list[tuple[str, str]] = []
    missing_mw_kind: list[str] = []

    for key, data in lib.items():
        smi = str(data.get("smiles") or "").strip()
        if smi:
            if Chem.MolFromSmiles(smi) is None:
                smiles_failed.append((key, smi))
                if log:
                    logger.warning(
                        "부형제 %s: smiles 파싱 실패 — RDKit이 읽지 못했습니다: %r", key, smi
                    )
            else:
                smiles_ok += 1

        for probe in data.get("reactive_probes") or []:
            probe_s = str(probe)
            if Chem.MolFromSmiles(probe_s) is None:
                probe_failed.append((key, probe_s))
                if log:
                    logger.warning(
                        "부형제 %s: reactive_probe 파싱 실패: %r", key, probe_s
                    )

        kind = data.get("mw_kind")
        if kind is None:
            missing_mw_kind.append(key)
            if log:
                logger.warning("부형제 %s: mw_kind 미지정 — 몰비 계산이 방어되지 않습니다", key)
        elif kind not in MW_KIND_VALUES:
            bad_mw_kind.append((key, str(kind)))
            if log:
                logger.warning(
                    "부형제 %s: mw_kind 값이 잘못됐습니다 %r (허용: %s)",
                    key, kind, ", ".join(MW_KIND_VALUES),
                )

    return {
        "total": len(lib),
        "smiles_ok": smiles_ok,
        "smiles_failed": smiles_failed,
        "probe_failed": probe_failed,
        "bad_mw_kind": bad_mw_kind,
        "missing_mw_kind": missing_mw_kind,
    }


def mw_kind_for(excipient_key: str, library: dict[str, dict[str, Any]] | None = None) -> str:
    """부형제 키 → `mw_kind`. 미등재·미지정이면 "undefined" (안전한 쪽).

    p17이 몰비를 계산해도 되는지 판단할 때 씁니다.
    """
    lib = COMPLETE_EXCIPIENT_MASTER if library is None else library
    key = LEGACY_EXCIPIENT_ALIASES.get(excipient_key, excipient_key)
    kind = (lib.get(key) or {}).get("mw_kind")
    return kind if kind in MW_KIND_VALUES else MW_KIND_UNDEFINED


# 레거시 키 → 신규 키
#
# ⚠️ 2026-09-03 — 이 표가 **부형제 별칭의 단일 출처(single source of truth)**
#    입니다. p72._NITRITE_MATRIX_ALIASES 가 같은 매핑 4건을 따로 들고 있어
#    한쪽만 고치면 어긋났습니다. 이제 p72가 이 표를 import 해 병합합니다.
#    새 별칭은 **여기에만** 추가하세요.
LEGACY_EXCIPIENT_ALIASES: dict[str, str] = {
    "Lactose": "Lactose_Monohydrate",
    "PEG_Polyethylene_Glycol": "Polyethylene_Glycol_4000",
    "PVP_Crospovidone": "Crospovidone",
    "Polyethylene_Glycol_PEG": "Polyethylene_Glycol_4000",
    "Crospovidone_PVP": "Crospovidone",
    "Microcrystalline_Cellulose_MCC": "Microcrystalline_Cellulose",
    "Sodium_Lauryl_Sulfate_SLS": "Sodium_Lauryl_Sulfate",
}
