"""화면 다국어 — English (기본) · 한국어.

ICH M7 관련 용어는 임의로 번역하지 않았습니다. 전부 가이드라인 원문
표현이며 근거는 `md파일/ICH_M7_공식용어집.md`에 있습니다.

사용법 (템플릿):
    {{ t('nav.single') }}
    {{ t('single.scan_button') }}

사용법 (파이썬):
    from ..i18n import t
    flash(t('flash.draw_structure'), 'warn')

번역이 없는 키는 키 문자열을 그대로 돌려주고 로그에 남깁니다 —
누락을 조용히 넘기지 않기 위함입니다.
"""

from __future__ import annotations

import logging
import re
from typing import Any

from flask import g, has_request_context, session

log = logging.getLogger("prime.i18n")

LANGUAGES = {"en": "English", "ko": "한국어"}
DEFAULT_LANG = "en"
SESSION_KEY = "prime_lang"

_missing: set[str] = set()


# ---------------------------------------------------------------------------
# 번역 사전
# ---------------------------------------------------------------------------

STRINGS: dict[str, dict[str, str]] = {
    # ===================== 공통 =====================
    "app.name": {"en": "PRIME-DS", "ko": "PRIME-DS"},
    "app.tagline": {
        "en": "Predictive Risk &amp; Impurity Modeling<br>Engine for Drug Substance",
        "ko": "원료의약품 불순물 예측 및<br>ICH M7 위해 평가 플랫폼",
    },
    "app.guideline": {
        "en": "ICH M7(R2) · RDKit · Ketcher",
        "ko": "ICH M7(R2) · RDKit · Ketcher",
    },
    "app.disclaimer": {
        "en": "For research and screening purposes. Expert review is required for "
              "GMP and regulatory decisions.",
        "ko": "연구·스크리닝 참고용입니다. GMP·규제 최종 판단은 전문가 검토가 필요합니다.",
    },

    "common.yes": {"en": "Yes", "ko": "예"},
    "common.no": {"en": "No", "ko": "아니오"},
    "common.none": {"en": "None", "ko": "없음"},
    "common.save": {"en": "Save", "ko": "저장"},
    "common.run": {"en": "Run", "ko": "실행"},
    "common.clear": {"en": "Clear", "ko": "지우기"},
    "common.clear_results": {"en": "Clear results", "ko": "결과 지우기"},
    "common.back": {"en": "Back", "ko": "돌아가기"},
    "common.open": {"en": "Open", "ko": "열기"},
    "common.selected": {"en": "Selected", "ko": "선택됨"},
    "common.search": {"en": "Search", "ko": "검색"},
    "common.reset": {"en": "Reset", "ko": "초기화"},
    "common.apply": {"en": "Apply", "ko": "적용"},
    "common.download": {"en": "Download", "ko": "다운로드"},
    "common.preview": {"en": "Preview", "ko": "미리보기"},
    "common.cases": {"en": "", "ko": "건"},
    "common.count_suffix": {"en": "", "ko": "건"},
    "common.optional": {"en": "optional", "ko": "선택"},
    "common.structure": {"en": "Structure", "ko": "구조"},
    "common.smiles": {"en": "SMILES", "ko": "SMILES"},
    "common.mw": {"en": "MW", "ko": "분자량"},
    "common.formula": {"en": "Formula", "ko": "분자식"},
    "common.mechanism": {"en": "Mechanism", "ko": "메커니즘"},
    "common.pathway": {"en": "Pathway", "ko": "경로"},
    "common.category": {"en": "Category", "ko": "카테고리"},
    "common.name": {"en": "Name", "ko": "이름"},
    "common.loading": {"en": "Working…", "ko": "연산 중…"},
    "common.loading_sub": {"en": "This may take a moment.", "ko": "잠시만 기다려 주세요."},
    "common.copied": {"en": "Copied", "ko": "복사했습니다"},

    # ===================== 사이드바 =====================
    "nav.section.assessment": {"en": "Assessment", "ko": "예측"},
    "nav.section.search": {"en": "Search", "ko": "검색"},
    "nav.section.tools": {"en": "Tools", "ko": "도구"},

    "nav.single": {"en": "Single Compound", "ko": "단일 화합물"},
    "nav.process": {"en": "Process Impurities", "ko": "공정 불순물"},
    "nav.degradation": {"en": "Degradation Products", "ko": "분해산물"},
    "nav.excipient": {"en": "Excipient Interaction", "ko": "부형제 상호작용"},
    "nav.drug_search": {"en": "Drug Database", "ko": "의약품 DB 검색"},
    "nav.structure_search": {"en": "Structure Search", "ko": "구조 검색"},
    "nav.reports": {"en": "Reports", "ko": "보고서 보관함"},
    # 2026-09-30 수정: 값에 HTML 엔티티가 박혀 있어 Jinja 자동 이스케이프를 거치며
    #   화면에 "Settings &amp;amp; Index" 로 나왔습니다. 평문 & 로 되돌립니다.
    "nav.settings": {"en": "Settings & Index", "ko": "설정 · 인덱스"},
    "nav.about": {"en": "Engine Info", "ko": "엔진 정보"},

    "nav.stat.indexed": {"en": "Drug index", "ko": "의약품 색인"},
    "nav.stat.empirical": {"en": "Empirical DB", "ko": "실측 독성 DB"},
    "nav.stat.excipients": {"en": "Excipients", "ko": "부형제"},

    # ===================== 단일 화합물 =====================
    "single.title": {"en": "Single Compound Assessment", "ko": "단일 화합물 스캔"},
    "single.subtitle": {
        "en": "ICH M7(R2) hazard assessment and physicochemical properties for one structure",
        "ko": "구조 하나에 대한 ICH M7 위해 평가와 물성 예측",
    },
    "single.input_title": {"en": "Structure Input", "ko": "구조 입력"},
    "single.input_hint": {
        "en": "Draw the structure, or paste from ChemDraw with Ctrl+V",
        "ko": "그리거나, ChemDraw에서 복사해 캔버스에 Ctrl+V",
    },
    "single.target": {"en": "Query structure", "ko": "평가 대상 구조"},
    "single.target_tag": {"en": "Target", "ko": "평가 대상"},
    "single.reference": {"en": "Reference drug substance", "ko": "대조 기준 구조"},
    "single.reference_tag": {
        "en": "ICH M7 Class 4 · optional",
        "ko": "ICH M7 Class 4 · 선택",
    },
    "single.scan_button": {"en": "Run assessment", "ko": "독성 스캔 시작"},
    "single.busy": {"en": "Running assessment…", "ko": "독성 스캔 중…"},
    "single.busy_sub": {
        "en": "Empirical database lookup, structural alert scan, and statistical model.",
        "ko": "Layer A 실측 조회 · Layer B 구조 경보 · Layer C 통계 모델을 순서대로 실행합니다.",
    },
    "single.slow_scan": {"en": "Allow exhaustive database scan", "ko": "느린 전수 검색 허용"},
    "single.slow_scan_hint": {
        "en": "If the local index misses, query PubChem and scan the source CSV files. "
              "May take tens of seconds.",
        "ko": "자체 색인에서 못 찾으면 PubChem을 거쳐 원본 CSV까지 뒤집니다. 수십 초 걸릴 수 있습니다.",
    },
    "single.reference_note": {
        "en": "The reference structure is used only for <b>Class 4</b> assignment — when the "
              "same structural alert is present in the drug substance and that material has "
              "been tested and found non-mutagenic. It may be left empty.",
        "ko": "대조 기준 구조는 <b>Class 4 판정</b>에만 쓰입니다 — 평가 대상과 같은 구조 경보가 "
              "원료의약품에도 있고, 그 원료가 변이원성 음성으로 확인된 경우 등급을 완화합니다. "
              "비워 두어도 됩니다.",
    },
    "single.empty_title": {"en": "No assessment yet", "ko": "아직 스캔하지 않았습니다"},
    "single.empty_text": {
        "en": "Draw a structure above and select «Run assessment».",
        "ko": "위 캔버스에 구조를 그리고 «독성 스캔 시작»을 누르세요.",
    },

    # 판정 요약
    "result.summary": {"en": "Assessment Summary", "ko": "판정 요약"},
    "result.classification": {
        "en": "ICH M7(R2) classification",
        "ko": "ICH M7(R2) 유전독성 등급",
    },
    "result.escalated": {
        "en": "Escalated by statistical model",
        "ko": "통계 층 격상",
    },
    "result.preliminary": {"en": "Preliminary", "ko": "예비 판정"},
    "result.acceptable_intake": {
        "en": "Recommended acceptable intake",
        "ko": "권고 관리 기준",
    },
    "result.control_strategy": {"en": "Control strategy", "ko": "관리 전략"},
    "result.empirical_match": {"en": "Empirical database match", "ko": "실측 DB 조회"},
    "result.match_found": {"en": "Match", "ko": "일치"},
    "result.no_record": {"en": "No record", "ko": "기록 없음"},
    "result.ames_evidence": {
        "en": "Bacterial mutagenicity evidence",
        "ko": "세균 변이원성 근거",
    },
    "result.carc_evidence": {"en": "Carcinogenicity evidence", "ko": "발암성 근거"},
    "result.model_probability": {
        "en": "P(mutagenic) — statistical model",
        "ko": "통계 층 P(변이원성)",
    },
    "result.rationale": {"en": "Rationale", "ko": "판정 근거"},
    "result.decision_log": {"en": "Assessment log", "ko": "판정 근거 로그"},
    "result.lines": {"en": "lines", "ko": "줄"},
    "result.classification_failed": {
        "en": "Classification failed",
        "ko": "ICH M7 판정 실패",
    },
    "result.classification_failed_text": {
        "en": "The structure is valid but the assessment engine returned no result.",
        "ko": "구조는 유효하지만 판정 엔진이 결과를 내지 못했습니다.",
    },
    "result.heavy_atoms": {"en": "heavy atoms", "ko": "중원자"},

    # 상태 배지
    "badge.empirically_detected": {
        "en": "Empirically detected",
        "ko": "실측 검출 이력 있음",
    },
    "badge.theoretical": {"en": "Theoretical prediction", "ko": "이론적 예측"},
    "badge.alert_detected": {"en": "detected", "ko": "검출"},
    "badge.no_alert": {"en": "None detected", "ko": "미검출"},
    "badge.validated_short": {"en": "Detected", "ko": "실측"},
    "badge.theoretical_short": {"en": "Predicted", "ko": "이론"},

    # 구조 경보
    "alerts.title": {
        "en": "Expert rule-based structural alerts",
        "ko": "Layer B — Ashby-Tennant 구조 경보",
    },
    "alerts.subtitle": {"en": "Ashby-Tennant alert set", "ko": "Ashby-Tennant 알러트"},
    "alerts.none_title": {"en": "No structural alerts", "ko": "구조 경보 없음"},
    # 2026-10-06 수정: 15 → 23. 이 카드가 쓰는 p62.scan_ashby_tennant_alerts 는
    #   ASHBY_TENNANT_ALERTS 21종 + 알고리즘 PAH 2종 = 23종을 훑는다.
    #   p65 의 METABOLIZATION_FRAGMENT_SMARTS 4종은 앙상블 판정에는 들어가지만
    #   이 카드에는 나오지 않으므로 여기 숫자에 더하지 말 것 (About 는 27종).
    "alerts.none_text": {
        "en": "None of the 23 Ashby-Tennant / Benigni-Bossa patterns matched.",
        "ko": "Ashby-Tennant · Benigni-Bossa 경보 23종 중 일치하는 것이 없습니다.",
    },
    # 2026-09-30 수정: 하드코딩된 특이도가 52.7% 로 남아 있었는데 이 프로젝트의
    #   어떤 측정값과도 일치하지 않습니다 (검증셋 69.0 · 학습분할 59.9 · 랜덤분할 64.9).
    #   출처: 논문자료/02_검증결과/external_tools/scaffold_comparison_report.txt
    #   규칙 계층을 다시 측정하면 이 숫자도 함께 고쳐야 합니다.
    "alerts.highlight_note": {
        "en": "Atoms matching the alert are <span class='hl-red'>highlighted in red</span> in "
              "the structure image. A structural alert alone does not establish mutagenicity — "
              "the rule layer alone has a specificity of 69.0% on the scaffold-separated "
              "validation set (n = 1,288).",
        "ko": "구조식 그림에서 <span class='hl-red'>붉게 표시된 부분</span>이 경보에 해당하는 "
              "원자입니다. 경보가 걸렸다고 반드시 변이원성인 것은 아닙니다 — 규칙 계층 단독 특이도는 "
              "스캐폴드 분할 검증셋(n = 1,288)에서 69.0% 입니다.",
    },

    # 물성
    "props.title": {"en": "Physicochemical Properties", "ko": "물리화학적 물성"},
    "props.subtitle": {
        "en": "RDKit descriptors and in-house estimates",
        "ko": "RDKit 표준 계산 + 자체 추정",
    },
    "props.stage1": {
        "en": "Stage 1 — RDKit descriptors (reliable)",
        "ko": "1단계 — RDKit 표준 (신뢰 가능)",
    },
    "props.stage2": {"en": "Stage 2 — in-house estimates", "ko": "2단계 — 자체 알고리즘"},
    "props.metric": {"en": "Property", "ko": "지표"},
    "props.value": {"en": "Value", "ko": "값"},
    "props.stability": {"en": "Process stability", "ko": "공정 안정성"},
    "props.estimate_warning_title": {
        "en": "pKa, LogD and boiling point are rule-assigned estimates",
        "ko": "pKa · LogD · 끓는점은 추정값입니다",
    },
    "props.estimate_warning_text": {
        "en": "A fixed value is assigned by functional group: 4.5 if an acidic group is present, "
              "9.5 if a basic nitrogen is present, otherwise 7.0. These are not computed from "
              "structure and <b>must not be cited in publications or regulatory documents</b>. "
              "They are not used in the hazard assessment.",
        "ko": "산성기가 있으면 4.5, 염기성 질소가 있으면 9.5, 둘 다 없으면 7.0을 대입하는 방식입니다. "
              "구조로부터 계산한 값이 아니므로 <b>논문·규제 문서에 인용하지 마세요</b>. "
              "이 값들은 독성 판정에 전혀 사용되지 않습니다.",
    },

    # 3중 대조
    "crossmatch.title": {
        "en": "Cross-match — in-house database, PubChem, structural alerts",
        "ko": "3중 대조 — 자체 DB · PubChem · 구조 경보",
    },
    "crossmatch.inhouse": {"en": "In-house index", "ko": "자체 색인"},
    "crossmatch.opensource": {
        "en": "Open source (PubChem + Ashby-Tennant)",
        "ko": "오픈소스 (PubChem + Ashby-Tennant)",
    },
    "crossmatch.empty": {"en": "No cross-match results", "ko": "대조 결과 없음"},

    # Read-Across
    "readacross.title": {"en": "Read-Across — structural analogues", "ko": "Read-Across — 유사 구조"},
    "readacross.similarity": {"en": "Similarity", "ko": "유사도"},
    "readacross.source": {"en": "Source", "ko": "출처"},
    "readacross.unnamed": {"en": "Unnamed", "ko": "이름 없음"},
    "readacross.note": {
        "en": "Analogue test results are <b>supporting information only</b> and do not affect "
              "the assigned class.",
        "ko": "유사 구조의 시험 결과는 <b>참고 정보</b>입니다. 판정 등급에는 반영되지 않습니다.",
    },

    # OECD
    "oecd.title": {"en": "OECD QSAR Toolbox (AMBIT2)", "ko": "OECD QSAR (AMBIT2)"},
    "oecd.query": {"en": "Query", "ko": "조회"},
    "oecd.requery": {"en": "Query again", "ko": "다시 조회"},
    "oecd.busy": {"en": "Querying OECD QSAR…", "ko": "OECD QSAR 조회 중…"},
    "oecd.item": {"en": "Endpoint", "ko": "항목"},
    "oecd.result": {"en": "Result", "ko": "결과"},
    "oecd.intro": {
        "en": "Queries hydrolysis, biodegradability and mutagenicity alerts. "
              "Runs only when requested because it is a network call.",
        "ko": "가수분해 · 생분해 · 변이원성 알러트를 조회합니다. 네트워크 요청이라 버튼을 눌렀을 때만 실행됩니다.",
    },
    "oecd.proxy_title": {
        "en": "Calculations are performed by RDKit",
        "ko": "실제 계산은 RDKit입니다",
    },
    "oecd.proxy_text": {
        "en": "Even when the AMBIT bridge is running, <code>run_ambit_bridge.bat</code> starts an "
              "internal emulation server (p55), not the OECD QSAR Toolbox. All results come from "
              "an RDKit proxy, so <b>«OECD Toolbox was used» must not be stated</b> in "
              "publications or regulatory documents.",
        "ko": "AMBIT 서버가 떠 있어도 <code>run_ambit_bridge.bat</code>이 띄우는 것은 진짜 OECD Toolbox가 "
              "아니라 내부 흉내내기 서버(p55)입니다. 어느 경로로 가든 계산은 RDKit 프록시이므로 "
              "논문·규제 문서에 <b>«OECD Toolbox 사용»이라고 기술하면 안 됩니다</b>.",
    },
    "oecd.endpoint": {"en": "Endpoint", "ko": "엔드포인트"},
    "oecd.failed": {"en": "Query failed", "ko": "조회 실패"},

    # ===================== 공정 =====================
    "process.title": {"en": "Process Impurity Assessment", "ko": "공정 불순물 예측"},
    "process.subtitle": {
        "en": "Impurities arising from starting materials, reagents and process conditions",
        "ko": "출발물질·시약·공정 조건에서 생기는 유연물질",
    },
    "process.input_title": {"en": "Reaction Structures", "ko": "반응 구조 입력"},
    "process.input_hint": {
        "en": "At least one starting material is required",
        "ko": "출발물질은 최소 하나 필요합니다",
    },
    "process.run": {"en": "Run process assessment", "ko": "하이브리드 공정 예측 가동"},
    "process.busy": {"en": "Running process assessment…", "ko": "공정 예측 중…"},
    "process.busy_sub": {
        "en": "Reagent side reactions, NDSRI formation, degradation, PubChem and structural alerts.",
        "ko": "시약 부반응 · NDSRI · 분해산물 · PubChem · 구조 경보를 순서대로 확인합니다.",
    },
    "process.reactant1": {"en": "[1] Starting material 1", "ko": "[1] 출발물질 1"},
    "process.reactant2": {"en": "[2] Starting material 2", "ko": "[2] 출발물질 2"},
    "process.reagent": {"en": "[3] Reagent", "ko": "[3] 반응 시약"},
    "process.product": {"en": "[4] Drug substance", "ko": "[4] 최종 목적물"},
    "process.reference_toggle": {
        "en": "Reference drug substance (optional) — for ICH M7 Class 4",
        "ko": "대조 기준 구조 (선택) — ICH M7 Class 4 판정용",
    },
    "process.conditions": {"en": "Process Conditions", "ko": "공정 조건"},
    "process.master_reagents": {"en": "master reagents", "ko": "마스터 시약"},
    "process.stress": {"en": "Stress degradation pathway", "ko": "가혹 분해 경로"},
    "process.name_reaction": {"en": "Name reaction", "ko": "Name Reaction"},
    "process.name_reaction_hint": {
        "en": "Main and side products are simulated separately using rdchiral",
        "ko": "rdchiral로 주반응·부반응을 나눠 시뮬레이션합니다",
    },
    "process.solvents": {"en": "Solvents", "ko": "용매"},
    "process.acids": {"en": "Acids", "ko": "산"},
    "process.bases": {"en": "Bases", "ko": "염기"},
    "process.catalysts": {"en": "Catalysts &amp; initiators", "ko": "촉매·개시제"},
    "process.multi_hint": {"en": "Hold Ctrl to select multiple", "ko": "Ctrl 누르고 여러 개 선택"},
    "process.temperature": {"en": "Reaction temperature (°C)", "ko": "반응 온도 (°C)"},
    "process.time": {"en": "Reaction time (h)", "ko": "반응 시간 (h)"},
    "process.nitrite": {"en": "Nitrite", "ko": "아질산염"},
    "process.nitrite_label": {
        "en": "Nitrite present in process",
        "ko": "공정에 아질산염 포함",
    },
    "process.nitrite_hint": {
        "en": "Drawing a nitrite source in the reagent canvas activates this automatically",
        "ko": "시약 캔버스에 아질산염을 그리면 이 옵션 없이도 자동 활성화됩니다",
    },
    "process.nitrite_on": {"en": "Nitrosation pathway active", "ko": "아질산염 경로 ON"},

    "process.dual_warning_title": {
        "en": "This screen runs two independent assessment paths",
        "ko": "이 화면에는 계산 경로가 둘 있습니다",
    },
    "process.dual_warning_text": {
        "en": "Path A (hybrid engine) uses classifier p66 — <b>escalation only</b>, it never lowers "
              "a class. Path B (process impurity engine) uses classifier p67, which <b>may lower</b> "
              "a class. The same compound can therefore appear with different classes in the two "
              "tables. Published performance figures are based on Path A.",
        "ko": "경로 A(하이브리드 엔진)는 판정기 p66 — <b>등급을 올리기만</b> 합니다. "
              "경로 B(공정 불순물 엔진)는 판정기 p67 — <b>낮추기도</b> 합니다. "
              "같은 화합물이 두 표에서 다른 등급으로 나올 수 있습니다. 논문 수치는 경로 A 기준입니다.",
    },
    "process.path_a": {"en": "Path A — Hybrid Process Engine", "ko": "Ⓐ 하이브리드 공정 예측"},
    "process.path_a_sub": {
        "en": "p39 · classifier p66 (escalation only) — reference for published figures",
        "ko": "p39 · 판정 p66 (상향 전용) — 논문 기준",
    },
    "process.path_b": {"en": "Path B — Process Impurity Engine", "ko": "Ⓑ 공정 불순물 엔진"},
    "process.path_b_sub": {
        "en": "p33 · p34 · p35 · p36 · classifier p67 (may lower class)",
        "ko": "p33·p34·p35·p36 · 판정 p67 (하향 허용)",
    },
    "process.target_profile": {"en": "Drug Substance Profile", "ko": "최종 목적물 프로필"},
    "process.vulnerable_patterns": {"en": "Vulnerable patterns", "ko": "취약 패턴"},
    "process.rotatable": {"en": "Rotatable bonds", "ko": "회전결합"},
    "process.col_type": {"en": "Type", "ko": "구분"},
    "process.col_verdict": {"en": "Status", "ko": "판정"},
    "process.col_route": {"en": "Reaction / route", "ko": "반응/경로"},
    "process.col_alert": {"en": "Structural alert", "ko": "구조 경보"},
    "process.no_impurities": {"en": "No impurities predicted", "ko": "예측된 유연물질이 없습니다"},
    "process.predicted": {"en": "Predicted impurities", "ko": "예측 유연물질"},
    "process.ai_generated": {"en": "AI generated", "ko": "AI 생성"},
    "process.db_verified": {"en": "DB match ≥85%", "ko": "DB 검증 ≥85%"},
    "process.high_risk": {"en": "ICH M7 Class 1–3", "ko": "ICH M7 Class 1–3"},
    "process.path_b_empty_title": {
        "en": "Path B returned no records",
        "ko": "경로 B 결과가 0건입니다",
    },
    "process.path_b_empty_text": {
        "en": "With no name reaction and no reagents selected, this engine has no rules to apply. "
              "Select a name reaction such as Suzuki, or choose solvents, acids, bases or catalysts, "
              "then run again. Path A results are unaffected.",
        "ko": "Name Reaction을 «선택안함»으로 두고 시약도 고르지 않으면 이 엔진은 적용할 규칙이 없습니다. "
              "Suzuki 같은 Name Reaction을 고르거나 용매·산·염기·촉매를 선택한 뒤 다시 실행해 보세요. "
              "경로 A 결과는 위에 그대로 있습니다.",
    },
    "process.bridge_connected": {"en": "External AI bridge connected", "ko": "외부 AI 브릿지 연결됨"},
    "process.bridge_unavailable": {
        "en": "ASKCOS not connected — RDKit engine used",
        "ko": "ASKCOS 미연결 — RDKit 엔진으로 예측",
    },
    "process.degradation_title": {"en": "Process Degradation Products", "ko": "공정 분해산물"},
    "process.degradation_parent": {"en": "Degradation parent", "ko": "분해 대상 구조"},
    "process.parent_differs": {
        "en": "differs from entered drug substance",
        "ko": "입력한 목적물과 다름",
    },
    "process.parent_differs_hint": {
        "en": "The engine re-resolved the target from the prediction results",
        "ko": "엔진이 예측 결과로부터 대상을 다시 정했습니다",
    },
    "process.degradation_skipped_title": {
        "en": "Degradation products were not assessed",
        "ko": "분해산물은 예측하지 않았습니다",
    },
    "process.degradation_skipped_text": {
        "en": "This step is skipped when the stress degradation pathway is set to «none».",
        "ko": "가혹 분해 경로를 «선택안함»으로 두면 이 단계를 건너뜁니다.",
    },
    "process.xlsx_disabled": {
        "en": "Cannot be generated — no records",
        "ko": "결과가 0건이라 만들 수 없습니다",
    },

    # ===================== 분해산물 =====================
    "degradation.title": {"en": "Degradation Product Assessment", "ko": "분해산물 예측"},
    "degradation.subtitle": {
        "en": "Degradation products under stress conditions and multi-generation degradation tree",
        "ko": "가혹 조건에서 생기는 분해산물과 다세대 분해 트리",
    },
    "degradation.conditions": {"en": "Stress Conditions", "ko": "분해 조건"},
    "degradation.parent": {"en": "Drug substance", "ko": "원료의약품 구조"},
    "degradation.run": {"en": "Predict degradation products", "ko": "분해산물 예측"},
    "degradation.busy": {"en": "Predicting degradation products…", "ko": "분해산물 예측 중…"},
    "degradation.busy_sub": {
        "en": "Applying stress-specific reaction rules and scoring by database similarity.",
        "ko": "가혹 조건 규칙을 적용하고 DB 유사도로 신뢰도를 매깁니다.",
    },
    "degradation.stress_env": {"en": "Stress condition", "ko": "가혹 조건"},
    "degradation.stress_hint": {
        "en": "Different rule sets apply to each condition",
        "ko": "조건마다 적용되는 반응 규칙이 다릅니다",
    },
    "degradation.temperature": {"en": "Test temperature (°C)", "ko": "시험 온도 (°C)"},
    "degradation.temperature_hint": {
        "en": "Above 60 °C adds to the similarity score",
        "ko": "60°C를 넘으면 신뢰도에 가산됩니다",
    },
    "degradation.rules_title": {"en": "Rules applied", "ko": "적용되는 규칙"},
    "degradation.rules_text": {
        "en": "Hydrolytic conditions → ester, amide, lactam and carbamate cleavage · "
              "Oxidative conditions → N-oxide, S-oxide, epoxidation and 4 others · "
              "Photolytic → Norrish Type I and photo-oxidation",
        "ko": "가수분해 조건 → 에스터·아마이드·락탐·카바메이트 절단 / "
              "산화 조건 → N-oxide, S-oxide, 에폭시화 등 7종 / "
              "광분해 → Norrish Type I + 광산화",
    },
    "degradation.results": {"en": "Predicted Degradation Products", "ko": "예측 결과"},
    "degradation.count": {"en": "Degradation products", "ko": "예측 분해산물"},
    "degradation.parent_mw": {"en": "Parent MW", "ko": "Parent 분자량"},
    "degradation.avg_similarity": {"en": "Mean DB similarity", "ko": "평균 DB 매칭률"},
    "degradation.col_route": {"en": "Degradation route", "ko": "분해 경로"},
    "degradation.col_product": {"en": "Predicted product (SMILES)", "ko": "예측 분해산물 (SMILES)"},
    "degradation.col_similarity": {"en": "DB similarity", "ko": "신뢰도"},
    "degradation.similarity_note": {
        "en": "«DB similarity» is the Tanimoto similarity to the in-house database of ~13,000 "
              "structures. It indicates how closely the predicted structure resembles a known "
              "compound — <b>it is not experimental verification</b>.",
        "ko": "신뢰도는 자체 DB 13,000종과의 <b>구조 유사도</b>입니다. 실측 검증 결과가 아니라 "
              "\"이 구조가 알려진 화합물과 얼마나 닮았는가\"를 나타냅니다.",
    },
    "degradation.structures": {"en": "Degradation Product Structures", "ko": "분해산물 구조"},
    "degradation.empty_title": {"en": "No degradation products", "ko": "분해산물이 없습니다"},
    "degradation.empty_text": {
        "en": "No degradation route applies to this structure under the selected condition.",
        "ko": "이 구조와 조건 조합에서는 유의미한 분해 경로가 적용되지 않았습니다.",
    },
    "degradation.empty_text_process": {
        "en": "This structure has no reactive site for the selected stress condition.",
        "ko": "이 구조에는 선택한 가혹 조건이 적용되는 반응 자리가 없습니다.",
    },

    # 분해 트리
    "tree.title": {"en": "Multi-generation Degradation Tree", "ko": "다세대 분해 트리"},
    "tree.subtitle": {
        "en": "Tracks degradation products that degrade further",
        "ko": "분해산물이 또 분해되는 경로를 추적합니다",
    },
    "tree.build": {"en": "Build tree", "ko": "트리 생성"},
    "tree.busy": {"en": "Building degradation tree…", "ko": "분해 트리 생성 중…"},
    "tree.busy_sub": {
        "en": "Re-applying reaction rules generation by generation.",
        "ko": "세대를 거듭하며 반응 규칙을 재적용합니다.",
    },
    "tree.depth": {"en": "Maximum generations", "ko": "최대 세대"},
    "tree.depth_hint": {"en": "Generation 0 = drug substance", "ko": "Generation 0 = 원료의약품"},
    "tree.include_ai": {"en": "Statistical model", "ko": "AI 통합"},
    "tree.include_ai_label": {
        "en": "Show P(mutagenic) and ICH M7 class on each node",
        "ko": "각 노드에 P(변이원성) · ICH M7 등급 표시",
    },
    "tree.nodes": {"en": "nodes", "ko": "노드"},
    "tree.max_gen": {"en": "max", "ko": "최대"},
    "tree.generations": {"en": "generations", "ko": "세대"},
    "tree.graph_note": {
        "en": "Nodes can be dragged. Colour indicates ICH M7 class — red Class 1–2, "
              "amber Class 3, blue Class 4, green Class 5.",
        "ko": "노드를 끌어 옮길 수 있습니다. 색은 ICH M7 등급 — 빨강 Class 1·2, 노랑 3, 파랑 4, 초록 5.",
    },
    "tree.node_list": {"en": "Node list", "ko": "노드 목록"},
    "tree.no_nodes": {"en": "No nodes to display", "ko": "그릴 노드가 없습니다"},
    "tree.failed": {"en": "Tree generation failed", "ko": "트리 생성 실패"},

    # Backend Core API
    "api.title": {
        "en": "Backend Core API — properties and read-across as JSON",
        "ko": "Backend Core API — 물성 + Read-Across가 붙은 JSON",
    },
    "api.ph": {"en": "Formulation pH", "ko": "제형 pH"},
    "api.depth": {"en": "Maximum generations", "ko": "최대 세대"},
    "api.image": {"en": "Structure image", "ko": "구조 이미지"},
    "api.include_b64": {"en": "Include structure_image_base64", "ko": "structure_image_base64 포함"},
    "api.generate": {"en": "Generate API response", "ko": "API 응답 생성"},
    "api.busy": {"en": "Generating API response…", "ko": "API 응답 생성 중…"},
    "api.failed": {"en": "API response failed", "ko": "API 응답 실패"},

    # ===================== 부형제 =====================
    "excipient.title": {"en": "Excipient Interaction Assessment", "ko": "부형제 상호작용"},
    "excipient.subtitle": {
        "en": "Impurities formed between excipients and the drug substance in the drug product",
        "ko": "완제의약품에서 부형제와 원료가 만드는 불순물",
    },
    "excipient.input_title": {
        "en": "Drug Substance and Excipients",
        "ko": "원료의약품 + 부형제",
    },
    "excipient.api": {"en": "Drug substance", "ko": "원료의약품 (API)"},
    "excipient.catalog": {"en": "selectable", "ko": "선택 가능"},
    # 2026-09-30 신설: 선택 가능한 항목 수(143)와 상세 레코드 보유 수(76)는 서로 다른
    #   집합이라 화면에서 구분해 표시합니다.
    "excipient.records": {"en": "detailed records", "ko": "상세 레코드"},
    "excipient.master": {"en": "master", "ko": "마스터"},
    "excipient.opensource": {"en": "open source", "ko": "오픈소스"},
    "excipient.run": {"en": "Assess interactions", "ko": "상호작용 예측"},
    "excipient.busy": {"en": "Assessing interactions…", "ko": "상호작용 예측 중…"},
    "excipient.busy_sub": {
        "en": "Checking Maillard, formylation and nitrosation pathways.",
        "ko": "Maillard · 포르밀화 · 니트로소화 경로를 확인합니다.",
    },
    "excipient.selected": {"en": "Selected excipients", "ko": "선택한 부형제"},
    "excipient.select_prompt": {
        "en": "Choose from the list below.",
        "ko": "아래 목록에서 골라 주세요.",
    },
    "excipient.scope_title": {
        "en": "Reaction pathways implemented",
        "ko": "적용되는 반응 경로",
    },
    "excipient.scope_text": {
        "en": "Maillard imine formation (aldehyde + amine) · formylation (carboxylic acid + amine) · "
              "nitrosation, plus acylation, Michael addition, bisulfite adduct and peroxide "
              "oxidation. Dissolution kinetics are not simulated.",
        "ko": "Maillard 이민 (알데하이드 + 아민) · 포르밀화 (카복실산 + 아민) · "
              "니트로소화 · 아실화 · Michael 부가 · 중아황산 부가 · 과산화물 산화. 용출 속도론은 미구현입니다.",
    },
    "excipient.selection_title": {"en": "Excipient Selection", "ko": "부형제 선택"},
    "excipient.clear_selection": {"en": "Clear selection", "ko": "선택 초기화"},
    "excipient.name_search": {"en": "Search by name", "ko": "이름 검색"},
    "excipient.add": {"en": "Add to selection", "ko": "선택 목록에 추가"},
    "excipient.showing": {"en": "shown of", "ko": "종 표시 중 · 전체"},
    "excipient.select_all": {"en": "Select all shown", "ko": "표시된 항목 전체 선택"},
    "excipient.expand_all": {"en": "Expand all", "ko": "모두 펼치기"},
    "excipient.collapse_all": {"en": "Collapse all", "ko": "모두 접기"},
    "excipient.name": {"en": "Excipient", "ko": "부형제"},
    "excipient.filtered": {"en": "matching the filter", "ko": "필터 일치"},
    "excipient.of_total": {"en": "of", "ko": "/ 전체"},
    "excipient.no_match": {
        "en": "No excipient matches the filter",
        "ko": "조건에 맞는 부형제가 없습니다",
    },
    "excipient.caution": {"en": "Caution", "ko": "주의"},
    "excipient.nitrite_title": {
        "en": "Nitrite residue profiles — excipients eligible for the nitrosation pathway",
        "ko": "아질산염 잔류 프로파일 — 니트로소화 경로를 타는 부형제",
    },
    "excipient.nitrite_note": {
        "en": "Excipients not listed here <b>do not enter the nitrosation pathway</b>; only "
              "Maillard and formylation apply. Residue levels are stored as text ranges and are "
              "not used in any numerical calculation.",
        "ko": "여기 없는 부형제는 <b>니트로소화 경로를 타지 않습니다</b>. Maillard·포르밀화만 적용됩니다. "
              "잔류 수준은 문자열 구간이며 수치 계산에 쓰이지 않습니다.",
    },
    "excipient.results": {"en": "Interaction Results", "ko": "상호작용 결과"},
    "excipient.total": {"en": "Total candidates", "ko": "전체 후보"},
    "excipient.amine_sites": {"en": "Amine sites in drug substance", "ko": "API 아민 자리"},
    "excipient.secondary": {"en": "secondary", "ko": "2급"},
    "excipient.tertiary": {"en": "tertiary", "ko": "3급"},
    "excipient.nitrosatable": {"en": "Nitrosatable", "ko": "니트로소화 가능"},
    "excipient.not_nitrosatable": {"en": "Not nitrosatable", "ko": "니트로소화 불가"},
    "excipient.rule": {"en": "Rule", "ko": "규칙"},
    "excipient.nitrite_exposure": {"en": "Nitrite", "ko": "아질산염"},
    "excipient.alpha_h_present": {"en": "α-hydrogen present", "ko": "α-수소 있음"},
    "excipient.alpha_h_absent": {"en": "No α-hydrogen — mitigated", "ko": "α-수소 없음 → 완화"},
    "excipient.alpha_h_note": {
        "en": "A nitrosamine without an α-hydrogen cannot undergo metabolic activation and is "
              "reduced to Class 5 by the chemical validity layer.",
        "ko": "α-수소가 없는 니트로사민은 대사활성화가 어려워 Class 5로 완화됩니다 (p74 화학적 타당성 검사).",
    },
    "excipient.hint_class": {
        "en": "Preliminary class estimate from the excipient engine",
        "ko": "부형제 엔진의 예비 등급 추정입니다",
    },
    # --- 1:1 이원 배합적합성 ---
    "binary.mode": {"en": "Study design", "ko": "평가 방식"},
    "binary.mode_binary": {
        "en": "Binary 1:1 (w/w) — one excipient at a time",
        "ko": "이원 1:1 (w/w) — 부형제 한 쌍씩",
    },
    "binary.mode_matrix": {
        "en": "Full formulation — all excipients together",
        "ko": "전체 배합 — 모든 부형제 동시",
    },
    "binary.mode_hint": {
        "en": "Binary 1:1 (w/w) mixtures are the standard preformulation compatibility "
              "design. Evaluating one pair at a time also prevents an excipient from being "
              "masked when another produces the same impurity.",
        "ko": "이원 1:1 (w/w) 혼합물은 전처방 배합적합성 시험의 표준 설계입니다. "
              "한 쌍씩 보면 같은 불순물을 만드는 다른 부형제에 가려지는 문제도 없어집니다.",
    },
    "binary.ratio": {"en": "Mass ratio (drug : excipient)", "ko": "배합비 (원료 : 부형제)"},
    "binary.ratio_hint": {
        "en": "1:1 by weight is the conventional binary compatibility ratio",
        "ko": "중량 기준 1:1이 이원 배합적합성의 관례적 비율",
    },
    "binary.results": {"en": "Binary Compatibility Results", "ko": "이원 배합적합성 결과"},
    "binary.ranking": {"en": "Excipient Risk Ranking", "ko": "부형제 위험도 순위"},
    "binary.ranking_sub": {
        "en": "Highest-concern pair first",
        "ko": "가장 주의가 필요한 조합부터",
    },
    "binary.rank": {"en": "Rank", "ko": "순위"},
    "binary.score": {"en": "Risk score", "ko": "위험 점수"},
    "binary.band": {"en": "Band", "ko": "구간"},
    "binary.band_high": {"en": "High", "ko": "높음"},
    "binary.band_moderate": {"en": "Moderate", "ko": "중간"},
    "binary.band_low": {"en": "Low", "ko": "낮음"},
    "binary.band_none": {"en": "No interaction", "ko": "상호작용 없음"},
    "binary.mole_ratio": {"en": "Molar ratio", "ko": "몰비"},
    "binary.mole_ratio_hint": {
        "en": "Moles of excipient per mole of drug substance at the stated mass ratio",
        "ko": "위 중량비에서 원료 1몰당 부형제 몰 수",
    },
    "binary.species": {"en": "Species", "ko": "생성물"},
    "binary.mechanisms": {"en": "Mechanisms", "ko": "반응 기전"},
    "binary.pairs": {"en": "Pairs evaluated", "ko": "평가한 쌍"},
    "binary.interacting": {"en": "With interaction", "ko": "상호작용 있음"},
    "binary.clean": {"en": "No interaction", "ko": "이상 없음"},
    "binary.total_species": {"en": "Total species", "ko": "전체 생성물"},
    "binary.heatmap": {"en": "Interaction Heat Map", "ko": "상호작용 히트맵"},
    "binary.heatmap_sub": {
        "en": "Excipient × reaction mechanism · cell value = number of species",
        "ko": "부형제 × 반응 기전 · 숫자는 생성물 건수",
    },
    "binary.clean_title": {"en": "No Interaction Detected", "ko": "상호작용이 확인되지 않은 조합"},
    "binary.clean_note": {
        "en": "These pairs were evaluated and no reaction rule applied. This is an "
              "<b>assessed result</b>, not an untested combination — but the rule set covers "
              "only Maillard, formylation and nitrosation, so it is not a general "
              "statement of compatibility.",
        "ko": "평가했으나 적용되는 반응 규칙이 없었던 조합입니다. <b>평가한 결과</b>이지 "
              "미시험이 아닙니다 — 다만 규칙은 Maillard·포르밀화·니트로소화 세 갈래뿐이므로 "
              "배합적합성 전반에 대한 보증은 아닙니다.",
    },
    "binary.detail": {"en": "Pair detail", "ko": "쌍 상세"},
    "binary.collapsed": {
        "en": "positional isomers collapsed",
        "ko": "위치 이성질체 접힘",
    },
    "binary.collapsed_note": {
        "en": "Where one reaction rule produced many products differing only in which "
              "hydroxyl or amine reacted, representative structures are shown and the "
              "remainder are counted. The transformation is chemically the same.",
        "ko": "한 반응 규칙이 어느 수산기·아민이 반응했는지만 다른 산물을 여러 개 만들 때, "
              "대표 구조만 보여주고 나머지는 건수로 접었습니다. 화학적으로 같은 변환입니다.",
    },
    "binary.rules": {"en": "Reaction rules applied", "ko": "적용 반응 규칙"},
    "binary.no_mw": {"en": "not defined", "ko": "정의 없음"},
    "binary.polymer_note": {
        "en": "Polymers and inorganic excipients have no single molecular weight, so a molar "
              "ratio is not defined. Reactivity is still evaluated through the reactive "
              "fragment recorded for each excipient.",
        "ko": "고분자·무기 부형제는 단일 분자량이 없어 몰비가 정의되지 않습니다. "
              "반응성은 부형제마다 기록된 반응 조각(reactive probe)으로 평가됩니다.",
    },
    "binary.score_note": {
        "en": "The risk score (0–100) prioritises which pairs to examine first. It is derived "
              "from the most severe species in the pair — impurity class, ICH M7 class hint, "
              "nitrite exposure — with a small bonus for additional species, and a reduction "
              "when no α-hydrogen is present. <b>It is not a regulatory limit.</b>",
        "ko": "위험 점수(0–100)는 어느 조합을 먼저 검토할지 정하는 우선순위 지표입니다. "
              "쌍에서 가장 심각한 생성물 한 건(불순물 분류 · ICH M7 등급 추정 · 아질산염 노출)을 "
              "기준으로 하고, 생성물이 여럿이면 소폭 가산, α-수소가 없으면 감산합니다. "
              "<b>규제 기준값이 아닙니다.</b>",
    },
    "binary.dedup_fix": {
        "en": "Why this differs from the full-formulation mode",
        "ko": "전체 배합 모드와 결과가 다른 이유",
    },
    "binary.dedup_fix_text": {
        "en": "In full-formulation mode all excipients are processed in one pass with a shared "
              "duplicate filter. When several excipients form the <b>same</b> impurity — as "
              "reducing sugars do — only the first one is reported and the others appear to "
              "have no interaction. Binary mode evaluates each pair independently, so every "
              "excipient is credited.",
        "ko": "전체 배합 모드는 모든 부형제를 한 번에 처리하면서 중복 필터를 공유합니다. "
              "여러 부형제가 <b>같은</b> 불순물을 만들면 — 환원당이 그렇습니다 — 먼저 처리된 "
              "하나만 보고되고 나머지는 상호작용이 없는 것처럼 보입니다. "
              "이원 모드는 쌍마다 독립적으로 평가하므로 전부 잡힙니다.",
    },

    "excipient.empty_title": {"en": "No interaction candidates", "ko": "상호작용 후보 없음"},
    "excipient.empty_text": {
        "en": "No reaction rule applies between the selected excipients and this structure.",
        "ko": "선택한 부형제와 이 구조 사이에 적용되는 반응 규칙이 없습니다.",
    },

    # ===================== 검색 =====================
    "search.drug_title": {"en": "Drug Database Search", "ko": "의약품 DB 검색"},
    "search.drug_subtitle": {
        "en": "Find a drug by name or CID and assess its related compounds",
        "ko": "이름 또는 CID로 찾고 관련 화합물을 평가합니다",
    },
    "search.drug_heading": {"en": "Drug Search", "ko": "의약품 검색"},
    "search.drug_hint": {
        "en": "Korean names work too — an alias table of 70 entries maps them to English",
        "ko": "한글 이름도 됩니다 — 별칭표 70개가 영문으로 바꿔 줍니다",
    },
    "search.drug_label": {"en": "Drug name or CID", "ko": "의약품명 · CID"},
    "search.results": {"en": "Search Results", "ko": "검색 결과"},
    "search.col_drug": {"en": "Drug", "ko": "의약품"},
    "search.col_file": {"en": "File", "ko": "파일"},
    "search.col_score": {"en": "Match", "ko": "일치도"},
    "search.not_found": {"en": "No results", "ko": "찾지 못했습니다"},
    "search.not_found_hangul": {
        "en": "The Korean name may not be in the alias table. Try the English name.",
        "ko": "한글로 검색하셨습니다. 별칭표에 없는 이름일 수 있으니 영문명으로 다시 시도해 보세요.",
    },
    "search.not_found_text": {
        "en": "Check the spelling, or search by CID.",
        "ko": "철자를 확인하거나 CID로 검색해 보세요.",
    },
    "search.related": {"en": "Related Compounds", "ko": "관련 화합물"},
    "search.of_total": {"en": "of", "ko": "행 중"},
    "search.rows": {"en": "rows", "ko": ""},
    "search.problem": {"en": "Search problem", "ko": "검색 중 문제"},

    "search.structure_title": {"en": "Structure Search", "ko": "구조 검색"},
    "search.structure_subtitle": {
        "en": "Substructure containment search and similarity top 3",
        "ko": "부분구조 포함 검색과 유사도 Top 3",
    },
    "search.query_structure": {"en": "Query structure", "ko": "검색할 구조"},
    "search.mode": {"en": "Substructure search scope", "ko": "부분구조 검색 범위"},
    "search.mode_fast": {
        "en": "<b>Fast</b> — representative compound per drug (milliseconds)",
        "ko": "<b>고속</b> — 의약품별 대표 화합물만 (밀리초)",
    },
    "search.mode_deep": {
        "en": "<b>Deep</b> — every CSV row (several minutes)",
        "ko": "<b>심층</b> — 전체 CSV 행 (수 분 소요)",
    },
    "search.run_substructure": {"en": "Substructure search", "ko": "부분구조 검색"},
    "search.run_similarity": {"en": "Similarity top 3", "ko": "유사도 Top 3"},
    "search.busy_substructure": {"en": "Searching substructures…", "ko": "부분구조 검색 중…"},
    "search.busy_substructure_sub": {
        "en": "Scanning the RDKit SubstructLibrary with multiple threads.",
        "ko": "RDKit SubstructLibrary를 멀티스레드로 훑습니다.",
    },
    "search.busy_similarity": {"en": "Computing similarity…", "ko": "유사도 계산 중…"},
    "search.busy_similarity_sub": {
        "en": "Morgan fingerprint Tanimoto comparison.",
        "ko": "Morgan 지문 Tanimoto 연산.",
    },
    "search.substructure_results": {"en": "Substructure Search Results", "ko": "부분구조 검색 결과"},
    "search.mode_deep_short": {"en": "Deep", "ko": "심층"},
    "search.mode_fast_short": {"en": "Fast", "ko": "고속"},
    "search.substructure_empty": {
        "en": "No containing structure found",
        "ko": "포함하는 구조를 찾지 못했습니다",
    },
    "search.substructure_empty_text": {
        "en": "Fast mode only searches representative compounds. Try deep mode.",
        "ko": "고속 모드는 의약품별 대표 화합물만 봅니다. 심층 모드로 다시 시도해 보세요.",
    },
    "search.similarity_results": {"en": "Similarity — Top 3", "ko": "구조 유사도 Top 3"},
    "search.similarity_sub": {"en": "Morgan fingerprint · Tanimoto", "ko": "Morgan 지문 · Tanimoto"},
    "search.similarity_empty": {"en": "No similar structures", "ko": "유사한 구조가 없습니다"},
    "search.similarity_empty_text": {
        "en": "The representative compound fingerprint index may be empty. Build it in Settings.",
        "ko": "대표 화합물 지문 인덱스가 비어 있을 수 있습니다. 설정에서 구축하세요.",
    },
    "search.col_compound": {"en": "Compound", "ko": "화합물"},
    "search.history": {"en": "Recent searches", "ko": "최근 검색"},
    "search.load": {"en": "Load", "ko": "불러오기"},
    "search.clear_history": {"en": "Clear history", "ko": "이력 지우기"},

    # ===================== 보고서 =====================
    "reports.title": {"en": "Reports", "ko": "보고서 보관함"},
    "reports.subtitle": {
        "en": "Reports available from the current session",
        "ko": "현재 세션에서 만들 수 있는 보고서",
    },
    "reports.regulatory": {"en": "Regulatory Report", "ko": "규제 보고서"},
    "reports.format": {"en": "ICH M7(R2) format", "ko": "ICH M7(R2) 양식"},
    "reports.single": {"en": "Single Compound Assessment", "ko": "단일 화합물 평가"},
    "reports.process": {"en": "Process Impurities", "ko": "공정 불순물"},
    "reports.degradation": {"en": "Degradation Products", "ko": "분해산물"},
    "reports.excipient": {"en": "Excipient Interaction", "ko": "부형제 상호작용"},
    "reports.drugdb": {"en": "Drug Database Screening", "ko": "의약품 DB 검색"},
    "reports.preview": {"en": "Preview", "ko": "미리보기"},
    "reports.ready": {"en": "Available", "ko": "결과 있음"},
    "reports.not_ready": {"en": "Not run yet", "ko": "아직 실행 안 함"},
    "reports.run_first": {
        "en": "Run the assessment on the",
        "ko": "먼저",
    },
    "reports.run_first_suffix": {"en": "screen first.", "ko": "화면에서 예측을 실행하세요."},
    "reports.note": {
        "en": "Reports are generated <b>when you click</b>. Large result sets may take a few "
              "seconds. Results are lost when the session is reset or after 8 hours.",
        "ko": "보고서는 <b>누른 시점에 생성</b>됩니다. 결과가 크면 몇 초 걸릴 수 있습니다. "
              "세션을 초기화하거나 8시간이 지나면 결과가 사라집니다.",
    },
    "reports.body": {"en": "Report body", "ko": "보고서 본문"},
    "reports.preview_title": {"en": "Report Preview", "ko": "보고서 미리보기"},
    "reports.empty": {"en": "Report body is empty", "ko": "본문이 비어 있습니다"},
    "reports.structures": {"en": "with structures", "ko": "구조식"},

    # ===================== 설정 =====================
    "settings.title": {"en": "Settings &amp; Index", "ko": "설정 · 인덱스"},
    "settings.subtitle": {
        "en": "Display options and search index management",
        "ko": "표시 옵션과 검색 인덱스를 관리합니다",
    },
    "settings.display": {"en": "Display Settings", "ko": "표시 설정"},
    "settings.max_show": {"en": "Maximum items shown", "ko": "최대 표시 개수"},
    "settings.max_show_hint": {
        "en": "How many result cards to display",
        "ko": "검색 결과 카드를 몇 개까지 보여줄지",
    },
    "settings.cols": {"en": "Cards per row", "ko": "한 줄 카드 수"},
    "settings.img_size": {"en": "Structure image size", "ko": "구조식 이미지 크기"},
    "settings.show_smiles": {"en": "SMILES", "ko": "SMILES 코드"},
    "settings.show_smiles_label": {
        "en": "Show SMILES strings on cards",
        "ko": "카드에 SMILES 문자열 표시",
    },
    "settings.language": {"en": "Language", "ko": "언어"},
    "settings.index": {"en": "Search Index", "ko": "검색 인덱스"},
    "settings.index_sub": {
        "en": "Determines structure search speed",
        "ko": "구조 검색 속도를 좌우합니다",
    },
    "settings.drug_index": {"en": "Drug index", "ko": "의약품 색인"},
    "settings.primary_fp": {"en": "Representative fingerprints", "ko": "대표 화합물 지문"},
    "settings.exact_index": {"en": "Exact match index", "ko": "Exact 인덱스"},
    "settings.substruct_index": {"en": "Substructure index", "ko": "하부구조 인덱스"},
    "settings.ready": {"en": "Ready", "ko": "준비됨"},
    "settings.absent": {"en": "Absent", "ko": "없음"},
    "settings.exact_empty_title": {
        "en": "The exact match index is empty",
        "ko": "Exact 인덱스가 비어 있습니다",
    },
    "settings.exact_empty_text": {
        "en": "The compounds table has 0 rows, so exact match search only covers representative "
              "compounds. «Build exact index» below will populate it, but it takes a long time.",
        "ko": "compounds 테이블이 0행이라 정확일치 검색이 대표 화합물로만 동작합니다. "
              "아래 «Exact 인덱스 구축»으로 채울 수 있지만 시간이 오래 걸립니다.",
    },
    "settings.reindex": {"en": "Refresh drug index", "ko": "색인 새로고침"},
    "settings.reindex_note": {"en": "Re-reads CSV files in data/", "ko": "data/ 폴더의 CSV를 다시 읽습니다"},
    "settings.build_primary": {"en": "Build representative fingerprints", "ko": "대표 화합물 지문 구축"},
    "settings.build_primary_note": {"en": "Required for similarity search", "ko": "유사도 검색에 필요"},
    "settings.build_substruct": {"en": "Build substructure index", "ko": "하부구조 인덱스 구축"},
    "settings.build_substruct_note": {"en": "Required for substructure search", "ko": "부분구조 검색에 필요"},
    "settings.build_compound": {"en": "Build exact index", "ko": "Exact 인덱스 구축"},
    "settings.build_compound_note": {"en": "Takes a long time", "ko": "장시간 소요"},
    "settings.index_busy": {"en": "Building index…", "ko": "인덱스 작업 중…"},
    "settings.session": {"en": "Session", "ko": "세션"},
    "settings.session_note": {
        "en": "Drawn structures, assessment results and excipient selections are held in server "
              "memory. Resetting clears all of them.",
        "ko": "캔버스에 그린 구조, 스캔 결과, 선택한 부형제가 모두 서버 메모리에 남아 있습니다. "
              "초기화하면 전부 지워집니다.",
    },
    "settings.reset_session": {"en": "Reset session", "ko": "세션 초기화"},
    "settings.reset_confirm": {
        "en": "Clear the entire session?",
        "ko": "현재 세션을 모두 지울까요?",
    },

    # ===================== 엔진 정보 =====================
    "about.title": {"en": "Engine Information", "ko": "엔진 정보"},
    "about.subtitle": {"en": "Components and known limitations", "ko": "구성 요소와 알려진 제약"},
    "about.tiers": {
        "en": "Four-tier ICH M7 assessment",
        "ko": "ICH M7 판정 4층 구조",
    },
    "about.tier": {"en": "Tier", "ko": "층"},
    "about.principle": {"en": "Principle", "ko": "원리"},
    "about.implementation": {"en": "Implementation", "ko": "구현"},
    "about.status": {"en": "Status", "ko": "상태"},
    "about.layer_a_principle": {"en": "Empirical database lookup", "ko": "실측 기록 조회"},
    "about.layer_b_principle": {
        "en": "Expert rule-based (Ashby-Tennant)",
        "ko": "전문가 규칙 (Ashby-Tennant)",
    },
    "about.layer_c_principle": {"en": "Statistical-based", "ko": "통계 학습"},
    # 2026-10-06: 키 이름이 «15» 로 굳어 있어 바꿨다. 앙상블의 규칙 계층은
    #   SMARTS 21종 + 알고리즘 PAH 2종 + 대사단편 4종으로 돈다.
    "about.alerts_impl": {
        "en": "21 SMARTS + 2 algorithmic + 4 metabolic-fragment patterns",
        "ko": "SMARTS 21종 + 알고리즘 2종 + 대사단편 4종",
    },
    "about.layer_d_principle": {
        "en": "Read-across arbitration",
        "ko": "읽기교차 중재",
    },
    "about.layer_d_status": {
        "en": "On B/C disagreement",
        "ko": "B·C 상충 시 발동",
    },
    "about.model": {"en": "model", "ko": "모델"},
    "about.normal": {"en": "Operational", "ko": "정상"},
    "about.error": {"en": "Error", "ko": "오류"},
    "about.model_missing": {"en": "Model file missing", "ko": "모델 파일 없음"},
    "about.never_run": {"en": "Never executed", "ko": "실행되지 않음"},
    "about.env_override_title": {
        "en": "PRIME_LAYER_C_MODEL environment variable is set",
        "ko": "PRIME_LAYER_C_MODEL 환경변수가 설정되어 있습니다",
    },
    "about.env_override_text": {
        "en": "model is currently forced. Unset it unless you are running an ablation experiment.",
        "ko": "모델이 강제 적용 중입니다. ablation 실험이 끝났다면 해제하세요.",
    },
    "about.escalation_title": {
        "en": "The consensus rule is escalation-only",
        "ko": "합산 규칙은 «상향 전용»입니다",
    },
    "about.escalation_text": {
        "en": "No tier can lower an assigned class. This is a deliberately conservative design "
              "that errs toward caution, at the cost of specificity (ensemble 67.0%).",
        "ko": "어느 층도 등급을 낮추지 못합니다. 안전 쪽으로 틀리도록 설계된 보수적 방식이며, "
              "그 대가로 특이도가 낮습니다 (앙상블 67.0%).",
    },
    "about.performance": {"en": "Validation Performance", "ko": "검증 성능"},
    "about.validation_set": {
        "en": "External validation set n = 1,288 · 588 positive · 700 negative "
              "(Murcko scaffold-separated)",
        "ko": "외부 검증셋 n = 1,288 · 양성 588 · 음성 700 (Murcko 스캐폴드 분할)",
    },
    "about.method": {"en": "Method", "ko": "방법"},
    "about.sensitivity": {"en": "Sensitivity", "ko": "민감도"},
    "about.specificity": {"en": "Specificity", "ko": "특이도"},
    "about.accuracy": {"en": "Accuracy", "ko": "일치율"},
    "about.mcc": {"en": "Matthews correlation coefficient", "ko": "Matthews 계수"},
    "about.layer_b_alone": {"en": "Expert rule-based alone", "ko": "Layer B 단독"},
    "about.layer_c_alone": {"en": "Statistical model alone", "ko": "Layer C 단독"},
    "about.ensemble": {"en": "Full ensemble", "ko": "전체 앙상블"},
    "about.performance_note": {
        "en": "The validation set does not overlap with the training set or the empirical "
              "database. Of 20 compounds escalated by the statistical layer alone "
              "(P > 0.70) despite carrying no structural alert, 15 were genuinely "
              "Ames-positive (precision 75.0%); read-across arbitration escalated a further "
              "31, of which 22 were positive (71.0%).",
        "ko": "검증셋은 학습셋 및 Layer A 실측 DB와 중복이 없습니다. "
              "구조 경보 없이 통계 계층이 단독으로 격상시킨 20건(P > 0.70) 가운데 15건이 실제 "
              "Ames 양성이었습니다 (정밀도 75.0%). 읽기교차 중재가 추가로 31건을 격상시켰고 "
              "그중 22건이 양성입니다 (71.0%).",
    },
    "about.limitations": {"en": "Known Limitations", "ko": "알려진 제약"},
    "about.limitations_sub": {
        "en": "Must not be cited in publications or regulatory documents",
        "ko": "논문·규제 문서에 인용하면 안 되는 항목",
    },
    "about.lim_retro_title": {
        "en": "Retrosynthesis and process economics figures",
        "ko": "역합성 · 공정경제성 지표",
    },
    "about.lim_retro_text": {
        "en": "Yield, bulk price, PMI and green chemistry scores are placeholders derived by "
              "mapping an MD5 hash of the structure onto a fixed range. No chemistry is involved.",
        "ko": "수율·원료가·PMI·그린케미스트리 점수는 구조식의 MD5 해시를 범위에 대응시킨 자리표시자입니다. "
              "화학적 계산이 전혀 없습니다.",
    },
    "about.lim_props_title": {"en": "pKa, LogD and boiling point", "ko": "pKa · LogD · 끓는점"},
    "about.lim_props_text": {
        "en": "One of three fixed values is assigned — 4.5 acidic, 9.5 basic, 7.0 neutral. "
              "This is classification, not prediction.",
        "ko": "산성기 4.5 / 염기성 9.5 / 중성 7.0 세 값 중 하나를 대입합니다. 예측이 아니라 분류입니다.",
    },
    "about.lim_oecd_title": {"en": "OECD Toolbox integration", "ko": "OECD Toolbox 연동"},
    "about.lim_oecd_text": {
        "en": "Even with the AMBIT server running, all calculations are performed by RDKit. "
              "«OECD Toolbox was used» must not be stated.",
        "ko": "AMBIT 서버가 떠 있어도 실제 계산은 RDKit입니다. «OECD Toolbox 사용»이라고 기술하면 안 됩니다.",
    },
    "about.lim_askcos_title": {"en": "ASKCOS and AiZynthFinder", "ko": "ASKCOS · AiZynthFinder"},
    "about.lim_askcos_text": {
        "en": "Only the client code is present; the engines themselves are not included. "
              "They function only when an external server is running.",
        "ko": "요청을 보내는 코드만 있고 엔진 자체는 포함되어 있지 않습니다. 외부 서버 구동 시에만 동작합니다.",
    },
    "about.lim_option_title": {
        "en": "«Option 3/4/5» labelling in reports",
        "ko": "보고서의 «Option 3/4/5» 표기",
    },
    "about.lim_option_text": {
        "en": "In ICH M7(R2), Options 1–4 are control strategy approaches, not assessment "
              "methods, and Option 5 does not exist. The labelling in report Section 3.4 is "
              "incorrect and needs to be corrected.",
        "ko": "ICH M7(R2)에서 Option 1~4는 평가 방법이 아니라 관리 전략이며 Option 5는 존재하지 않습니다. "
              "보고서 3.4절 표기가 잘못되어 있어 수정이 필요합니다.",
    },
    "about.composition": {"en": "Composition", "ko": "구성"},
    "about.degradation_version": {"en": "Degradation engine version", "ko": "분해 엔진 버전"},
    "about.excipient_catalog": {"en": "Excipient catalogue", "ko": "부형제 카탈로그"},
    "about.nitrite_profiles": {"en": "Nitrite residue profiles", "ko": "아질산염 프로파일"},
    "about.nitrite_profiles_value": {
        "en": "7 excipients — all others do not enter the nitrosation pathway",
        "ko": "7종 — 나머지 부형제는 니트로소화 경로를 타지 않습니다",
    },
    "about.excipient_rules": {"en": "Excipient reaction rules", "ko": "부형제 반응 규칙"},
    "about.excipient_rules_value": {
        "en": "Maillard imine and formylation, plus nitrosation",
        "ko": "Maillard 이민 · 포르밀화 2종 + 니트로소화",
    },
    "about.editor": {"en": "Structure editor", "ko": "구조 편집기"},
    "about.editor_value": {
        "en": "Ketcher (externally hosted — internet required)",
        "ko": "Ketcher (외부 호스팅 — 인터넷 필요)",
    },
    "about.total": {"en": "total", "ko": "총"},
    "about.records": {"en": "records", "ko": "건"},

    # ===================== 편집기 =====================
    "editor.empty": {"en": "Empty", "ko": "비어 있음"},
    "editor.saved": {"en": "Saved", "ko": "저장됨"},
    "editor.placeholder": {
        "en": "Draw a structure — the SMILES will appear here",
        "ko": "구조를 그리면 SMILES가 여기에 표시됩니다",
    },
    "editor.clear": {"en": "Clear canvas", "ko": "캔버스 비우기"},
    "editor.no_preview": {"en": "No preview", "ko": "미리보기 없음"},
    "editor.pasting": {"en": "Applying paste…", "ko": "붙여넣기 반영 중…"},

    # ===================== flash 메시지 =====================
    "flash.draw_structure": {
        "en": "Please draw a structure first.",
        "ko": "먼저 구조를 그려 주세요.",
    },
    "flash.parse_failed": {
        "en": "The structure could not be interpreted.",
        "ko": "구조를 해석하지 못했습니다.",
    },
    "flash.parse_failed_named": {
        "en": "could not be interpreted.",
        "ko": "의 구조를 해석하지 못했습니다.",
    },
    "flash.results_cleared": {"en": "Results cleared.", "ko": "결과를 지웠습니다."},
    "flash.select_stress": {
        "en": "Select a stress condition to apply degradation routes.",
        "ko": "가혹 조건을 선택해야 분해 경로가 적용됩니다.",
    },
    "flash.run_degradation_first": {
        "en": "Run the degradation assessment first, or enter a structure.",
        "ko": "먼저 분해산물 예측을 실행하거나 구조를 입력하세요.",
    },
    "flash.structure_required": {"en": "A structure is required.", "ko": "구조가 필요합니다."},
    "flash.select_excipient": {
        "en": "Select at least one excipient.",
        "ko": "부형제를 하나 이상 선택하세요.",
    },
    "flash.no_excipient_hits": {
        "en": "No interaction candidates. No reaction rule applies between the selected "
              "excipients and this structure.",
        "ko": "상호작용 후보가 없습니다. 선택한 부형제와 이 구조 사이에 적용되는 반응 규칙이 없습니다.",
    },
    "flash.need_starting_material": {
        "en": "Draw at least one starting material.",
        "ko": "출발물질을 최소 하나는 그려 주세요.",
    },
    "flash.need_product": {
        "en": "A drug substance or reference structure is needed to generate a report.",
        "ko": "최종 목적물이나 대조 구조 중 하나는 있어야 보고서를 만들 수 있습니다.",
    },
    "flash.scan_first": {
        "en": "Run the assessment first.",
        "ko": "먼저 구조를 스캔하세요.",
    },
    "flash.display_saved": {"en": "Display settings saved.", "ko": "표시 설정을 저장했습니다."},
    "flash.session_reset": {
        "en": "Session reset. Canvases and previous results have been cleared.",
        "ko": "세션을 초기화했습니다. 캔버스와 이전 결과가 모두 지워졌습니다.",
    },
    "flash.history_cleared": {"en": "Search history cleared.", "ko": "검색 이력을 지웠습니다."},
    "flash.language_changed": {"en": "Language changed.", "ko": "언어를 변경했습니다."},
    "flash.index_refreshed": {"en": "Index refreshed —", "ko": "색인 갱신 —"},
    "flash.registered": {"en": "registered", "ko": "등록"},
    "flash.skipped": {"en": "skipped", "ko": "건너뜀"},
    "flash.failed": {"en": "failed", "ko": "실패"},
    "flash.index_failed": {"en": "Index refresh failed:", "ko": "색인 갱신 실패:"},
    "flash.fp_built": {"en": "Representative fingerprints —", "ko": "대표 화합물 지문 —"},
    "flash.fp_failed": {"en": "Fingerprint build failed:", "ko": "지문 구축 실패:"},
    "flash.substruct_built": {"en": "Substructure index —", "ko": "하부구조 인덱스 —"},
    "flash.substruct_error": {"en": "Substructure index error:", "ko": "하부구조 인덱스 오류:"},
    "flash.substruct_failed": {"en": "Substructure index failed:", "ko": "하부구조 인덱스 실패:"},
    "flash.exact_built": {"en": "Exact index —", "ko": "Exact 인덱스 —"},
    "flash.exact_failed": {"en": "Exact index failed:", "ko": "Exact 인덱스 실패:"},
    "flash.new": {"en": "new", "ko": "신규"},
    "flash.scanned_rows": {"en": "rows scanned", "ko": "스캔 행"},
    "flash.files": {"en": "files", "ko": "파일"},
    "flash.report_ctx_failed": {
        "en": "Could not build the report context — TXT only.",
        "ko": "보고서 컨텍스트를 만들지 못해 TXT만 제공합니다.",
    },
    "flash.report_failed": {"en": "generation failed:", "ko": "생성 실패:"},
    "flash.report_failed_txt": {
        "en": "— download the TXT version instead.",
        "ko": "— TXT로 대신 내려받으세요.",
    },
    "flash.no_records_excel": {
        "en": "Path B returned no records, so the Excel report cannot be generated. Select a name "
              "reaction or reagents and run again.",
        "ko": "공정 불순물 엔진(경로 B) 결과가 0건이라 Excel을 만들 수 없습니다. "
              "Name Reaction이나 시약을 선택한 뒤 다시 실행해 보세요.",
    },
    "flash.empty_excel": {"en": "The generated Excel file is empty.", "ko": "Excel 생성 결과가 비어 있습니다."},
    "flash.no_scan": {
        "en": "No assessment result. Run a single compound assessment first.",
        "ko": "스캔 결과가 없습니다. 먼저 단일 화합물 스캔을 실행하세요.",
    },
    "flash.no_degradation": {"en": "No degradation results.", "ko": "분해산물 결과가 없습니다."},
    "flash.no_tree": {"en": "No degradation tree results.", "ko": "분해 트리 결과가 없습니다."},
    "flash.no_excipient": {"en": "No excipient interaction results.", "ko": "부형제 상호작용 결과가 없습니다."},
    "flash.no_drugdb": {
        "en": "No drug database results. Search for a drug and open a file first.",
        "ko": "의약품 DB 검색 결과가 없습니다. 먼저 의약품을 검색하고 파일을 선택하세요.",
    },
    "flash.no_process": {"en": "No process assessment results.", "ko": "공정 예측 결과가 없습니다."},

    # JS 문구
    "js.draw_first": {
        "en": "Please draw a structure first —",
        "ko": "구조를 먼저 그려 주세요 —",
    },
    "js.copy_failed": {"en": "Copy failed", "ko": "복사에 실패했습니다"},
}


# ---------------------------------------------------------------------------
# API
# ---------------------------------------------------------------------------


def current_lang() -> str:
    """현재 언어. 요청 컨텍스트 밖에서는 기본값."""
    if not has_request_context():
        return DEFAULT_LANG
    cached = getattr(g, "_prime_lang", None)
    if cached:
        return cached
    lang = session.get(SESSION_KEY) or DEFAULT_LANG
    if lang not in LANGUAGES:
        lang = DEFAULT_LANG
    g._prime_lang = lang
    return lang


def set_lang(lang: str) -> str:
    if lang not in LANGUAGES:
        lang = DEFAULT_LANG
    session[SESSION_KEY] = lang
    session.permanent = True
    if has_request_context():
        g._prime_lang = lang
    return lang


def t(key: str, *, lang: str | None = None, **kwargs: Any) -> str:
    """번역 문자열. 없으면 키를 그대로 돌려주고 한 번만 경고합니다."""
    language = lang or current_lang()
    entry = STRINGS.get(key)
    if entry is None:
        if key not in _missing:
            _missing.add(key)
            log.warning("번역 키 누락: %s", key)
        return key
    # 2026-09-30 수정: `or` 연쇄가 빈 문자열을 falsy 로 보아, 일부러 비워 둔 번역
    #   (영문 common.count_suffix 등)에서 키 문자열이 화면에 그대로 나왔습니다.
    #   None 여부로만 판정합니다.
    text = entry.get(language)
    if text is None:
        text = entry.get(DEFAULT_LANG)
    if text is None:
        text = key
    if kwargs:
        try:
            text = text.format(**kwargs)
        except (KeyError, IndexError, ValueError):
            pass
    return text


# ---------------------------------------------------------------------------
# 엔진 상수의 화면 표시 라벨
#
# 가혹 조건 이름은 엔진(p43·p33)에 한글로 박혀 있고 값 자체가 식별자로
# 쓰입니다. 값은 그대로 두고 화면에 보일 때만 영문으로 바꿉니다.
# 용어는 ICH Q1A/Q1B 표준 표현을 따랐습니다.
# ---------------------------------------------------------------------------

STRESS_LABELS: dict[str, str] = {
    # p43 — 분해산물 화면
    "선택 안함": "None",
    "산성 가혹 가수분해 (Acidic Hydrolysis)": "Acidic hydrolysis",
    "염기성 가혹 가수분해 (Alkaline Hydrolysis)": "Alkaline hydrolysis",
    "수분 노출 (Neutral Hydrolysis)": "Neutral hydrolysis (humidity)",
    "과산화 산화 조건 (Oxidative Stress)": "Oxidative stress (peroxide)",
    "광분해 조건 (Photolytic Degradation)": "Photolytic degradation",
    # p33 — 공정 화면
    "선택안함 (No Accelerated Stress)": "None (no accelerated stress)",
    "Acidic Hydrolysis (Harsh)": "Acidic hydrolysis (harsh)",
    "Alkaline Hydrolysis (Harsh)": "Alkaline hydrolysis (harsh)",
    "Oxidative Stress Conditions": "Oxidative stress conditions",
    # p34 · p13 — 선택 안함 계열
    "선택안함 (No Name Reaction)": "None (no name reaction)",
    "선택안함 (Not Selected)": "Not selected",
    "표준 가수분해/산화": "Standard hydrolysis / oxidation",
}


# p50_physicochemical.predict_chemical_stability가 돌려주는 한글 문구
ENGINE_TEXT: dict[str, str] = {
    "✅ 안정 (Stable)": "✅ Stable",
    "ℹ️ 보통 (Moderately Stable)": "ℹ️ Moderately stable",
    "⚠️ 불안정 (Unstable)": "⚠️ Unstable",
    "예측 불가": "Not predicted",
    "물성 데이터 부족": "Insufficient property data",
    "안정성 연산 오류": "Stability calculation error",
    "주요 가수분해·산화 취약 구조가 발견되지 않았으며, 화학적 거동이 안정할 것으로 예측됨.":
        "No major hydrolytic or oxidative liability was identified; the compound is "
        "predicted to be chemically stable.",
    "수용액 공정 조건에서 친수성 에스테르 결합의 급격한 가수분해(Hydrolysis) 리스크가 높음.":
        "High risk of rapid hydrolysis of the hydrophilic ester bond under aqueous "
        "process conditions.",
    "강산/강염기 혹은 고온 조건에서 제한적인 아미드 가수분해 가능성 있음.":
        "Limited amide hydrolysis is possible under strong acid, strong base or "
        "elevated temperature.",
    "추후 업데이트 예정": "To be updated",

    # p39 — 하이브리드 공정 엔진의 행 구분·판정 라벨
    "교차반응 / NDSRI": "Cross-reaction / NDSRI",
    "시약 부반응 / NDSRI": "Reagent side reaction / NDSRI",
    "미반응 출발물질·시약": "Unreacted starting material / reagent",
    "목적물 분해산물": "Degradation product of drug substance",
    "📌 [실측 데이터 기반 유연물질 (Validated)]": "📌 Empirically detected impurity",
    "🧪 [오픈소스 규칙 기반 예측 불순물 (Theoretical)]": "🧪 Theoretical prediction (rule-based)",
    "변이원성 BioAssay 미검출": "No mutagenicity bioassay found",
    "변이원성 시험 매핑 있음": "Mutagenicity assay data available",
    "⚪ PubChem 조회 실패": "⚪ PubChem lookup failed",
    "분류 실패": "Classification failed",
    "없음": "None",

    # p40 — 분해 취약성 패턴
    "에스터 결합 — 가수분해 취약": "Ester bond — susceptible to hydrolysis",
    "아미드 결합 — 산/염기 가수분해 취약": "Amide bond — susceptible to acid/base hydrolysis",
    "락탐 고리 — 개환 가수분해": "Lactam ring — ring-opening hydrolysis",
    "β-락탐 — 고리 개환": "β-Lactam — ring opening",
    "카르바메이트 — 가수분해": "Carbamate — hydrolysis",
    "설폰아미드 — 가수분해": "Sulfonamide — hydrolysis",
    "이미드 — 가수분해": "Imide — hydrolysis",
    "헤미아세탈 — 산성 가수분해": "Hemiacetal — acid hydrolysis",
    "이황화 결합 — 환원 분해": "Disulfide bond — reductive cleavage",
    "니트로기 — 환원 취약": "Nitro group — susceptible to reduction",
    "방향족 아민 — N-산화 취약": "Aromatic amine — susceptible to N-oxidation",
    "페놀 — 산화 취약": "Phenol — susceptible to oxidation",
    "주요 가수분해/산화 취약 패턴 미검출":
        "No major hydrolytic or oxidative liability detected",
    "공정 안정성": "Process stability",

    # p32 · p30 — 반응 규칙 라벨
    "미반응 출발물질 1": "Unreacted starting material 1",
    "미반응 출발물질 2": "Unreacted starting material 2",
    "미반응 반응 시약": "Unreacted reagent",
    "출발물질 1": "Starting material 1",
    "출발물질 2": "Starting material 2",
    "목적물": "Drug substance",
    "에스터 가수분해": "Ester hydrolysis",
    "아미드 가수분해": "Amide hydrolysis",
    "락탐 가수분해": "Lactam hydrolysis",
    "카르바메이트 가수분해": "Carbamate hydrolysis",
    "설폰아미드 가수분해": "Sulfonamide hydrolysis",
    "이미드 가수분해": "Imide hydrolysis",
    "API 가수분해 분해산물": "Hydrolytic degradation product",
    "API 자동산화 분해산물": "Autoxidation degradation product",
    "API NDSRI 분해": "NDSRI formation from drug substance",
    "알데하이드 + 2급 아민 → 환원아민화 부산물":
        "Aldehyde + secondary amine → reductive amination by-product",
    "카르복실산 + 2급 아민 → 아미드": "Carboxylic acid + secondary amine → amide",
    "2급 아민 + 아질산염 → N-니트로소amine":
        "Secondary amine + nitrite → N-nitrosamine",
    "2급 아민 + 시약 → N-니트로소amine": "Secondary amine + reagent → N-nitrosamine",
    "3급 아민 + 시약 → N-니트로소amine": "Tertiary amine + reagent → N-nitrosamine",
    "2급 아민 + 아질산염 → N-니트로소amine (NDSRI)":
        "Secondary amine + nitrite → N-nitrosamine (NDSRI)",
    "2급 아민 알킬화 부반응 (알킬화 시약)":
        "Secondary amine alkylation side reaction (alkylating reagent)",
    "방향족 N-산화 (산화제 시약)": "Aromatic N-oxidation (oxidising reagent)",
    "에스터 가수분해 (산성 시약)": "Ester hydrolysis (acidic reagent)",
    "카르보닐 환원 부반응 (환원 시약)": "Carbonyl reduction side reaction (reducing reagent)",
    "아질산염 소스 (NDSRI 위험)": "Nitrite source (NDSRI risk)",
    "알킬화 시약": "Alkylating reagent",
    "산화제": "Oxidising agent",
    "환원제": "Reducing agent",
    "산성 시약": "Acidic reagent",
    "아질산": "Nitrite",

    # p68 — 3중 대조 표 (컬럼명 · 값)
    "항목": "Item",
    "결과": "Result",
    "매칭 결과": "Match result",
    "매칭 화합물명": "Matched compound",
    "화합물명": "Compound",
    "분자식 / MW": "Formula / MW",
    "의약품 파일": "Source file",
    "자체 DB (13,040종)": "In-house database (13,040 entries)",
    "자체 DB에 동일 Canonical SMILES 없음": "No matching canonical SMILES in the in-house database",
    "Ames/변이원성 BioAssay": "Bacterial mutagenicity bioassay",
    "Ashby-Tennant 스캔": "Ashby-Tennant alert scan",
    "검출 알러트": "Alerts detected",
    "PubChem 요약에 관련 시험 미검출": "No relevant assay found in the PubChem summary",
    "📌 판정: [실측 데이터 기반 유연물질 (Validated)]": "📌 Empirically detected impurity",
    "🧪 판정: [오픈소스 규칙 기반 예측 불순물 (Theoretical)]": "🧪 Theoretical prediction (rule-based)",
    "원료의약품 SMILES 해석 실패": "Could not parse the drug substance SMILES",

    # p42 · p43 — 분해 경로 · 산물 명칭
    "SMIRKS 에스테르 가수분해": "SMIRKS ester hydrolysis",
    "SMIRKS 아미드 가수분해": "SMIRKS amide hydrolysis",
    "SMIRKS 카바메이트 가수분해": "SMIRKS carbamate hydrolysis",
    "SMIRKS 3급 아민 N-Oxide": "SMIRKS tertiary amine N-oxide",
    "SMIRKS 티오에테르 자동산화 (Sulfoxide)": "SMIRKS thioether autoxidation (sulfoxide)",
    "SMIRKS 에테르 탈메틸화 (Autoxidation)": "SMIRKS ether demethylation (autoxidation)",
    "범용 에스테르 가수분해 (Universal Ester Cleavage)": "Universal ester cleavage",
    "범용 아미드 가수분해 (Universal Amide Cleavage)": "Universal amide cleavage",
    "락탐 고리 개환 가수분해 (Lactam Ring Opening)": "Lactam ring-opening hydrolysis",
    "카바메이트 가수분해 (Carbamate Hydrolysis)": "Carbamate hydrolysis",
    "N-자동산화 (N-Oxide Formation)": "N-autoxidation (N-oxide formation)",
    "S-자동산화 (S=O Formation)": "S-autoxidation (S=O formation)",
    "S-자동산화 (Sulfoxide Formation)": "S-autoxidation (sulfoxide formation)",
    "질소 자동산화 (N-Oxidation Fallback)": "Nitrogen autoxidation (N-oxidation)",
    "황 자동산화 (S-Oxidation Fallback)": "Sulfur autoxidation (S-oxidation)",
    "알데하이드 자동산화 (Aldehyde → Carboxylic acid)":
        "Aldehyde autoxidation (aldehyde → carboxylic acid)",
    "Norrish Type I 광화학 절단": "Norrish Type I photochemical cleavage",
    "광산화 (Photo-oxidation) —": "Photo-oxidation —",
    "광분해": "Photolysis",
    "분해 fragment": "Degradation fragment",
    "유기 분해 fragment": "Organic degradation fragment",
    "에스테르 절단 분해체": "Ester cleavage product",
    "아미드 절단 분해체": "Amide cleavage product",
    "카르복실산 산화 분해체": "Carboxylic acid oxidation product",
    "에폭사이드 산화 분해체": "Epoxide oxidation product",
    "설폭시드 산화체": "Sulfoxide oxidation product",
    "설폭시드형 산화체": "Sulfoxide-type oxidation product",
    "설폭시드/아울폭사이드형 산화체": "Sulfoxide / sulfone-type oxidation product",
    "N-옥사이드 자동산화체": "N-oxide autoxidation product",
    "N-옥사이드 대사체": "N-oxide metabolite",
    "N-옥사이드": "N-oxide",
    "N-탈알킬화 2차 아민 대사체": "N-dealkylated secondary amine metabolite",
    "다성분 하이드록실화 잔기": "Multi-site hydroxylation residue",
    "아릴아민 산화 전구체": "Arylamine oxidation precursor",
    "탈할로겐·탈수소 부산물 — 방향족 할로겐 제거체":
        "Dehalogenation / dehydrogenation by-product — aryl halide removal",
    "탈할로겐·탈수소 부산물": "Dehalogenation / dehydrogenation by-product",
    "에스테르 절단 — 카르복실산 잔기": "Ester cleavage — carboxylic acid residue",
    "카르복실산 잔기": "Carboxylic acid residue",
    "알코올/페놀릭 잔기": "Alcohol / phenolic residue",
    "아민/아마인 잔기": "Amine residue",
    "새로 생긴 방향족 OH": "Newly formed aromatic OH",
    "유사 화합물": "Analogous compound",
    "DB FP 인덱스 미매칭 — 메커니즘 연산만 적용":
        "No database fingerprint match — mechanism-based assessment only",
    "RDKit 구조 해석 실패 — DB 매칭 생략":
        "RDKit could not parse the structure — database matching skipped",
    "SMILES 없음 — DB 매칭 생략": "No SMILES — database matching skipped",
    "RAG(DB 유사도 조회, LLM 검색증강생성 아님)":
        "Database similarity lookup (not LLM retrieval-augmented generation)",

    # 신뢰도 등급 · 예상 비율 · 위험 신호
    "🔥 실제 검출 이력 일치 (매우 확실)": "🔥 Matches a documented detection (very high)",
    "🔥 실제 검출 이력과 구조 90% 이상 일치 (검증 완료)":
        "🔥 ≥90% structural match to a documented detection (verified)",
    "🔥 실제 검출 이력과 구조 90% 이상 일치 (신뢰도 높음)":
        "🔥 ≥90% structural match to a documented detection (high)",
    "⚠️ 이론적 메커니즘 일치 및 유사 작용기 DB 확인":
        "⚠️ Mechanism match with analogous functional groups in the database",
    "⚠️ 화학 메커니즘 및 유사 작용기 유효성 확인":
        "⚠️ Chemically plausible mechanism and functional group",
    "⚠️ 메커니즘적 유효": "⚠️ Mechanistically plausible",
    "💡 이론적 예측 최적 모델 기반 도출": "💡 Derived from the predictive model",
    "💡 정보화학 모델 기반 가상 시뮬레이션 결과": "💡 In silico simulation result",
    "💡 이론적 가능성": "💡 Theoretical possibility",
    "🔴 위험 (유전독성 경고 구조 포함)": "🔴 High concern (genotoxic alert present)",
    "🟡 주의 (잠재적 위험)": "🟡 Caution (potential concern)",
    "🟢 안전 (독성 위험 낮음)": "🟢 Low concern",
    "주요 (Major, ~40–55%)": "Major (~40–55%)",
    "부분 (Minor, ~20–35%)": "Minor (~20–35%)",
    "미량 (Trace, <15%)": "Trace (<15%)",
    "미량 (Trace, ~5–15%)": "Trace (~5–15%)",
    "미량 (Trace, ~8–15%)": "Trace (~8–15%)",
    "고온 가혹 (Accelerated)": "Accelerated (elevated temperature)",
    "중간 가혹 (Intermediate)": "Intermediate",
    "저온/실온 (Ambient)": "Ambient",

    # to_report_row 컬럼명
    "분해 경로": "Degradation route",
    "예측 분해산물 구조 (SMILES)": "Predicted degradation product (SMILES)",
    "실시간 연산 분자량 (g/mol)": "Molecular weight (g/mol)",
    "ICH M7 독성 등급": "ICH M7 class",
    "예측 신뢰도 및 DB 매칭률 (%)": "Database similarity (%)",
    "예측 분해산물명": "Predicted product name",
    "위험 신호": "Risk signal",
    "예상 비율": "Estimated proportion",
    "반응 메커니즘": "Reaction mechanism",
    "DB RAG 매칭": "Database match",
    "Ashby-Tennant 알러트": "Ashby-Tennant alerts",
}

# 문장 안에 섞여 나오는 조각. 긴 것부터 치환합니다.
ENGINE_FRAGMENTS: dict[str, str] = {
    # 아래 두 라벨은 2026-08-26에 폐기된 reagent_ndsri_user 규칙의 것입니다.
    # 이전 세션에서 만든 보고서에 남아 있을 수 있어 번역만 남겨 둡니다.
    "2급/3급 아민 + 시약 → NDSRI (요청 SMIRKS)":
        "Secondary/tertiary amine + reagent → NDSRI (retired rule)",
    "2급/1급 아민 + 아질산염 → N-니트로소amine (요청 SMIRKS)":
        "Primary/secondary amine + nitrite → N-nitrosamine (retired rule)",
    "2급 아민 + 시약(아질산) → N-니트로소amine":
        "Secondary amine + nitrite reagent → N-nitrosamine",
    "3급 아민 + 시약(아질산) → 니트로소 탈알킬화 → NDSRI":
        "Tertiary amine + nitrite reagent → nitrosative dealkylation → NDSRI",
    "3급 아민 + 시약 → 니트로소 탈알킬화 → NDSRI":
        "Tertiary amine + reagent → nitrosative dealkylation → NDSRI",
    "3급 아민 + 시약(아질산) → N-니트로소amine":
        "Tertiary amine + nitrite reagent → N-nitrosamine",
    "2급 아민 + 시약 → N-니트로소amine": "Secondary amine + reagent → N-nitrosamine",
    "3급 아민 + 시약 → N-니트로소amine": "Tertiary amine + reagent → N-nitrosamine",
    "니트로소 탈알킬화": "nitrosative dealkylation",
    "비변이원성 니트로사민 (ICH M7 Class 5 위험도 낮음)":
        "Non-mutagenic nitrosamine (ICH M7 Class 5 — low concern)",
    "비변이원성 니트로사민": "Non-mutagenic nitrosamine",
    "위험도 낮음": "low concern",
    "🔥 Positive (실측)": "🔥 Positive (empirical)",
    "Positive (실측)": "Positive (empirical)",
    "Negative (실측)": "Negative (empirical)",
    "Positive (위험)": "Positive",
    "Negative (안전)": "Negative",
    "실측": "empirical",
    "⚪ PubChem에 등록된 동일 구조 없음": "⚪ No matching structure in PubChem",
    "PubChem에 등록된 동일 구조 없음": "No matching structure in PubChem",
    "주요 가수분해/산화 취약 패턴 미검출":
        "No major hydrolytic or oxidative liability detected",
    "에스테르 가수분해": "Ester hydrolysis",
    "아미드 가수분해": "Amide hydrolysis",
    "락탐 개환": "Lactam ring opening",
    "주형 모핵 잔기": "Core scaffold residue",
    "소분자 파편": "Small-molecule fragment",
    "분해산물 잔기": "Degradation residue",
    "Exact SMILES 일치": "Exact SMILES match",
    "가혹 온도": "Stress temperature",
    "중간 가혹": "intermediate stress",
    "고온 가혹": "accelerated",
    "취약 패턴": "Liability patterns",
    "출발물질 1": "starting material 1",
    "출발물질 2": "starting material 2",
    "미반응": "unreacted",
    "다단계 분해": "multi-step degradation",
    "검출": "detected",
    "미검출": "not detected",
    "안전": "no concern",
    "위험": "concern",
    "없음": "none",
    "건": "",
}

# 접두어가 일치하면 그 뒤를 붙여 쓰는 동적 문구
ENGINE_PREFIX: tuple[tuple[str, str], ...] = (
    ("취약 패턴 ", "Liability patterns — "),
    ("다단계 분해 (Gen ", "Multi-step degradation (Gen "),
    ("PubChem Ames/변이원성 BioAssay ", "PubChem mutagenicity bioassays: "),
)

_FRAGMENTS_SORTED = sorted(ENGINE_FRAGMENTS.items(), key=lambda kv: -len(kv[0]))
_HANGUL = re.compile(r"[가-힣]")


def option_label(value: str, *, lang: str | None = None) -> str:
    """엔진 상수를 화면 표시용 문자열로. 값 자체는 바뀌지 않습니다."""
    language = lang or current_lang()
    if language != "en":
        return value
    text = str(value or "").strip()
    if not text:
        return text

    hit = STRESS_LABELS.get(text) or ENGINE_TEXT.get(text)
    if hit:
        return hit

    for prefix, replacement in ENGINE_PREFIX:
        if text.startswith(prefix):
            return _translate_fragments(replacement + text[len(prefix):])

    # "규칙명 (출발물질 1)" 처럼 조각이 섞인 문자열
    if _HANGUL.search(text):
        return _translate_fragments(text)

    return text



def _translate_fragments(text: str) -> str:
    """알려진 조각을 긴 것부터 바꿔 넣습니다."""
    out = text
    for ko, en in _FRAGMENTS_SORTED:
        if ko in out:
            out = out.replace(ko, en)
    for ko, en in ENGINE_TEXT.items():
        if ko in out:
            out = out.replace(ko, en)
    # 조각 제거로 생긴 이중 공백 정리
    return re.sub(r"\s{2,}", " ", out).strip()


def missing_keys() -> list[str]:
    """실행 중 발견된 누락 키 — 검증용."""
    return sorted(_missing)


def audit() -> dict[str, list[str]]:
    """사전 자체를 점검합니다. 번역이 빠진 항목을 찾습니다."""
    problems: dict[str, list[str]] = {"missing_en": [], "missing_ko": [], "identical": []}
    for key, entry in STRINGS.items():
        en = (entry.get("en") or "").strip()
        ko = (entry.get("ko") or "").strip()
        if key in INTENTIONALLY_EMPTY:
            continue
        if not en:
            problems["missing_en"].append(key)
        if not ko:
            problems["missing_ko"].append(key)
        if en and ko and en == ko and not _looks_intentionally_same(en):
            problems["identical"].append(key)
    return problems


# 한쪽이 비어 있는 것이 의도인 키
#   영어에는 «건» 같은 수량 접미사가 없고, 한국어에는 «rows» 같은 단위어가 없습니다.
INTENTIONALLY_EMPTY: frozenset[str] = frozenset({
    "common.cases",
    "common.count_suffix",
    "search.rows",
})


def _looks_intentionally_same(text: str) -> bool:
    """제품명·약어처럼 번역하지 않는 것이 맞는 경우."""
    return text in {
        "PRIME-DS",
        "SMILES",
        "Target",
        "ICH M7(R2) · RDKit · Ketcher",
        "ICH M7 Class 1–3",
    }
