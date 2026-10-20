"""
PRIME-DS 모듈 의존 지도 생성기.

코드베이스를 AST로 훑어 `PRIME-DS_모듈지도.html` 안의 DATA 블록을 다시 씁니다.
HTML의 화면 부분(CSS·JS)은 건드리지 않고 데이터만 교체하므로,
코드가 바뀌면 이 스크립트만 다시 돌리면 됩니다.

    python s10_build_module_map.py

모듈의 «역할» 한 줄과 «변경 이력»은 자동으로 뽑을 수 없으므로
아래 ROLES / HIST 에 손으로 적습니다. 새 모듈을 만들면 여기에 추가하세요.
(빠뜨리면 실행할 때 경고가 뜹니다.)
"""

from __future__ import annotations

import ast
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
HTML = os.path.join(ROOT, "PRIME-DS_모듈지도.html")

# ---------------------------------------------------------------------------
# 밴드 — 지도의 열. 왼쪽부터 파이프라인 순서.
# ---------------------------------------------------------------------------

BANDS = [
    {"id": "b0",  "tag": "0x",  "label": "앱 · 설정"},
    {"id": "bwv", "tag": "web", "label": "웹 화면 · 라우팅"},
    {"id": "bws", "tag": "web", "label": "웹 서비스 · 기반"},
    {"id": "b1",  "tag": "1x",  "label": "레지스트리 · 규칙"},
    {"id": "b2",  "tag": "2x",  "label": "데이터 · 검색"},
    {"id": "b3",  "tag": "3x",  "label": "반응 · 공정 예측"},
    {"id": "b4",  "tag": "4x",  "label": "분해 예측"},
    {"id": "b5",  "tag": "5x",  "label": "물성 · QSAR"},
    {"id": "b6",  "tag": "6x",  "label": "독성 · ICH M7"},
    {"id": "b7",  "tag": "7x",  "label": "부형제 · 니트로사민"},
    {"id": "b8",  "tag": "8x",  "label": "리포트"},
    {"id": "b9",  "tag": "9x",  "label": "산출 · 보조 앱"},
    {"id": "bs",  "tag": "s",   "label": "오프라인 스크립트"},
]

# ---------------------------------------------------------------------------
# 역할 — 카드에 표시되는 한 줄 설명
# ---------------------------------------------------------------------------

ROLES: dict[str, str] = {
    # 0x 앱 · 설정
    "p00_flask_app": "Flask 실행 진입점",
    "p01_prime_config": "로컬 설정 로더",
    "p02_degradation_constants": "캐시 버전 상수",
    "p03_branding": "브랜딩 · 로그 접두사",
    "p04_cache_paths": "캐시 경로",
    "p05_safeguards": "연산 세이프가드",
    # 1x 레지스트리 · 규칙
    "p13_chemical_registry": "시약 마스터 DB",
    "p14_dynamic_registry": "동적 시약 등록",
    "p15_process_pipeline": "공정 예측 파이프라인",
    "p16_ich_m7_terms": "ICH M7 공식 용어 상수",
    "p17_binary_compatibility": "1:1 이원 배합 평가",
    "p18_excipient_reaction_rules": "부형제 반응 규칙 14종",
    "p19_excipient_library_ext": "부형제 라이브러리 확장",
    # 2x 데이터 · 검색
    "p20_db_index": "파일 메타 인덱스",
    "p21_drug_search_aliases": "한글명 검색어 확장",
    "p22_drug_data": "의약품 CSV 로드",
    "p23_chem_search_engine": "구조 검색 엔진",
    "p24_db_smiles_lookup": "SMILES 실측 조회",
    "p25_pubchem_client": "PubChem 연동",
    "p26_import_index": "데이터 임포트",
    # 3x 반응 · 공정
    "p30_reaction_engine": "NDSRI SMIRKS 엔진",
    "p31_process_reaction_rules": "공정 반응 템플릿",
    "p32_process_reaction_engine": "2-반응물 예측",
    "p33_process_impurity_engine": "공정 유연물질 엔진",
    "p34_name_reaction_process": "Name Reaction 시뮬",
    "p35_universal_reaction_impurity_engine": "범용 부반응 엔진",
    "p36_suzuki_impurity_rules": "Suzuki 부반응 규칙",
    "p37_retro_synthesis_process": "역합성 · 공정경제",
    "p38_external_ai_reaction_bridge": "외부 AI 반응 브리지",
    "p39_hybrid_process_engine": "하이브리드 공정 예측",
    # 4x 분해
    "p40_degradation_engine": "취약 구조 스캔",
    "p41_degradation_smirks_rules": "분해 SMIRKS 카탈로그",
    "p42_degradation_rag_engine": "조건부 분해 RAG",
    "p43_degradation_predictor": "분해산물 예측",
    "p44_degradation_shelf_life_filter": "보관 안정성 여과",
    # 5x 물성 · QSAR
    "p50_physicochemical": "물성 계산",
    "p51_prime_backend_core": "백엔드 코어 스위트",
    "p52_ambit_client": "OECD QSAR 클라이언트",
    "p53_oecd_rdkit_fallback": "QSAR 로컬 대체",
    "p55_ambit_bridge_server": "AMBIT 브리지 서버",
    # 6x 독성 · ICH M7
    "p60_m7_empirical_db": "M7 실측 DB (Layer A)",
    "p61_build_m7_empirical_db": "실측 DB 빌더",
    "p62_toxicity_alerts": "구조 알러트 (Layer B)",
    "p63_deepchem_predictor": "Ames 통계 모델 (Layer C)",
    "p64_gnn_mutagenicity_predictor": "AttentiveFP GNN",
    "p65_m7_ensemble_matrix": "3층 앙상블 매트릭스",
    "p66_ich_m7": "ICH M7 합의 판정",
    "p67_ich_m7_classifier": "M7 Class 1~5 판정",
    "p68_hybrid_toxicity": "하이브리드 독성 평가",
    # 7x 부형제
    "p70_excipient_master_library": "부형제 마스터",
    "p71_excipient_opensource_registry": "FDA IID 레지스트리",
    "p72_excipient_interaction_engine": "부형제–API 상호작용",
    "p74_nitrosamine_chemical_validity": "니트로사민 타당성",
    # 8x 리포트
    "p80_formal_assessment_report": "규제용 평가 리포트",
    "p81_report_export": "DOCX · PDF 출력",
    "p82_pdf_typesetting": "PDF 조판 헬퍼",
    "p83_m7_executive_visuals": "M7 요약 시각화",
    "p84_structure_images": "구조식 이미지",
    "p86_degradation_report_export": "분해 리포트",
    "p87_degradation_tree_report_export": "분해 트리 PDF",
    "p88_excipient_report_export": "부형제 리포트",
    "p89_drug_db_report_export": "DB 검색 리포트",
    # 9x
    "p90_report_generator": "공정 리포트 생성",
    "p91_prime_integrated_report": "통합 JSON",
    "p93_weasyprint_support": "WeasyPrint 설정",
    "p96_read_across_flask_app": "Read-across API",
    "p97_build_attentivefp_model": "GNN 모델 빌더",
    # 스크립트
    "s1_setup_layer_c_ames": "Layer C 재학습",
    "s2_build_layer_a_from_ames": "Layer A DB 구축",
    "s3_fix_layer_a_carcinogens": "Layer A 오류 수정",
    "s4_recheck_validation_set": "검증셋 중복 점검",
    "s5_run_pipeline_validation": "파이프라인 외부검증",
    "s6_retrain_layer_c_model": "Layer C 경량 재학습",
    "s7_run_case_study": "니트로사민 케이스",
    "s8_run_excipient_case_study": "부형제 케이스",
    "s9_run_binary_compatibility_study": "이원 배합 케이스",
    "s10_build_module_map": "이 지도 생성기",
    # webapp
    "webapp/__init__": "Flask 앱 팩토리",
    "webapp/state": "서버측 세션",
    "webapp/filters": "템플릿 필터",
    "webapp/i18n": "화면 번역 사전 (EN/KO)",
    "webapp/excipient_terms": "부형제 영문명 · 주의문구",
    "webapp/services/__init__": "서비스 패키지",
    "webapp/services/chem": "SMILES · 물성 · M7 판정",
    "webapp/services/single": "단일 구조 스캔",
    "webapp/services/process": "공정 불순물 (경로 A·B)",
    "webapp/services/degradation": "분해산물 · 분해 트리",
    "webapp/services/excipient": "부형제 상호작용",
    "webapp/services/search": "의약품 · 구조 검색 · 이력",
    "webapp/views/__init__": "뷰 패키지",
    "webapp/views/main": "홈 · 설정 · 엔진 정보",
    "webapp/views/single": "단일 구조 화면",
    "webapp/views/process": "공정 불순물 화면",
    "webapp/views/degradation": "분해산물 화면",
    "webapp/views/excipient": "부형제 화면",
    "webapp/views/search": "검색 화면",
    "webapp/views/reports": "보고서 다운로드",
    "webapp/views/api": "구조 이미지 · 캔버스 API",
    "webapp/views/context": "템플릿 공통 컨텍스트",
}

# ---------------------------------------------------------------------------
# 변경 이력 — 상세 패널에 표시. 최신이 위로 간다.
#   kind : new(신규) · fix(수정) · note(참고) · open(미해결)
# ---------------------------------------------------------------------------

_FLASK = "2026-08-25"
_I18N = "2026-08-25"
_EXC = "2026-08-25"
_P30 = "2026-08-26"

HIST: dict[str, list[dict]] = {

    # ---------------------------------------------------- 0x
    "p00_flask_app": [
        {"d": _FLASK, "k": "new", "t": "신규 — Streamlit 진입점을 대체",
         "b": "이전에는 <code>p00_app.py</code>(81KB · 2,011줄)가 진입점이었다. "
              "실행 명령이 <code>streamlit run app.py</code>에서 "
              "<code>python p00_flask_app.py</code>로 바뀌었다. "
              "127.0.0.1:8501에서 브라우저가 자동으로 열리고, 부팅 시 "
              "RDKit · 모델 파일 · DB · PRIME_LAYER_C_MODEL 환경변수를 점검한다."},
        {"d": _FLASK, "k": "note", "t": "구버전 보관 위치",
         "b": "<code>_archive/streamlit_ui/</code> 에 p00_app · p10 · p11 · p12 · "
              "p54 · p85 · p92 · p95 와 ketcher_component 가 그대로 있다. "
              "되돌리려면 루트로 옮기고 <code>streamlit run p00_app.py</code>."},
        {"d": _FLASK, "k": "fix", "t": "구버전의 import 누락 버그를 넘겨받지 않음",
         "b": "<code>p00_app.py</code>는 <code>expand_search_queries</code>(L408)와 "
              "<code>has_hangul</code>(L1737)을 import 없이 쓰고 있었다. "
              "색인 미스 폴백과 검색 0건 분기에서 NameError가 나던 상태였다."},
    ],

    # ---------------------------------------------------- 1x
    "p15_process_pipeline": [
        {"d": _FLASK, "k": "new", "t": "신규 — p12에서 계산부만 뽑아냄",
         "b": "공정 불순물 예측의 핵심인 <code>run_process_impurity_prediction</code>이 "
              "화면 모듈 <code>p12_process_conditions_ui.py</code> 안에 있었다. "
              "이 파일은 최상단에서 <code>import streamlit</code>을 하므로 "
              "Streamlit 없이는 호출할 수 없었다. 로직은 그대로 두고 화면 의존성만 걷어냈다."},
        {"d": _FLASK, "k": "fix", "t": "시약 선택이 조용히 무시되던 버그",
         "b": "<code>build_process_reagent_context_from_ids</code>는 <b>소문자</b> 키만 받는다 "
              "(<code>solvents · acids · bases · catalysts</code>). 대문자로 넘기면 "
              "예외 없이 빈 결과를 돌려준다. 용매를 골라도 반영되지 않던 원인. "
              "<code>build_process_conditions</code>에서 소문자로 통일했다."},
        {"d": _FLASK, "k": "fix", "t": "시약 컨텍스트 키 이름 불일치",
         "b": "반환 dict의 키는 <code>primary_solvent · solvent_props</code>인데 "
              "<code>solvent · solvent_bp_c</code>로 읽고 있었다. "
              "용매 이름이 계속 «—»로 표시되던 원인."},
    ],
    "p16_ich_m7_terms": [
        {"d": _I18N, "k": "new", "t": "신규 — ICH M7 공식 용어 상수",
         "b": "Table 1의 Definition과 Proposed action을 <b>가이드라인 원문 그대로</b> 넣었다. "
              "화면·보고서에서 문자열을 직접 적으면 표현이 갈라지므로 반드시 이 모듈을 쓴다. "
              "근거는 <code>md파일/ICH_M7_공식용어집.md</code>."},
        {"d": _I18N, "k": "note", "t": "Control Options는 1~4번만 존재",
         "b": "§8.1의 Option 1~4는 <b>원료의약품 관리 전략</b>이지 평가 방법이 아니다. "
              "Option 5는 존재하지 않는다. TTC 기본값 1.5 µg/day."},
    ],
    "p17_binary_compatibility": [
        {"d": _EXC, "k": "new", "t": "신규 — 1:1 (w/w) 이원 배합 평가",
         "b": "전처방 배합적합성 시험 관례에 맞춰 부형제를 <b>한 쌍씩</b> 평가한다. "
              "부형제마다 엔진을 새로 만들어 <code>seen_smiles</code> 공유를 끊었다. "
              "위험도 = 가장 심각한 한 건 + 건수 가산점(최대 15), 100 상한. "
              "밴드 high ≥60 · moderate ≥30 · low >0."},
        {"d": _EXC, "k": "fix", "t": "전체 배합 방식의 위음성을 드러냄",
         "b": "기존 방식은 부형제 사이에 중복 제거 집합을 공유해 "
              "<b>같은 불순물을 만드는 부형제가 여럿이면 하나만 남았다.</b> "
              "부형제 76종 전수 확인 결과 원료 4종에서 37건 → 198건 (+161), "
              "원료마다 24~32종이 통째로 가려져 있었다."},
        {"d": "2026-08-26", "k": "open", "t": "논문 수치 갱신 필요",
         "b": "p30 수정 후 재측정했다. 부형제 143종 전수로 "
              "<b>634건 → 1,021건 (+387), 원료마다 48~49종이 가려진다.</b> "
              "다만 기존 문서의 대표 사례였던 «메트포민 + PEG 4000이 상호작용 없음» 은 "
              "<b>더 이상 재현되지 않는다</b> (전체배합 0→2건). 산물 구성이 바뀌면서 "
              "중복 충돌이 그 자리에서 일어나지 않게 됐다. "
              "표본이 작으면 사례가 흔들리므로 논문에는 5종 패널 대신 전수 수치를 쓰는 편이 낫다."},
    ],
    "p18_excipient_reaction_rules": [
        {"d": _EXC, "k": "new", "t": "신규 — 반응 규칙 2종 → 14종",
         "b": "Maillard(이민·엔아민) · 포르밀화 · 아마이드 축합 · 에스터화 · 에스터교환 · "
              "에스터 아미놀리시스 · Michael 부가(아민·티올) · 아황산부가물(중성·음이온) · "
              "과산화물 N-산화 · S-산화 · 헤미아세탈."},
        {"d": _EXC, "k": "fix", "t": "기존 Maillard 규칙이 화학적으로 틀렸음",
         "b": "<code>[NX3;H2,H1]</code> 하나로 처리해 <b>아마이드 질소에 이민을 붙였다</b> "
              "(아테놀올 → <code>CC(=O)N=CC(O)CO</code>). 2급 아민에서 중성 이민을 "
              "만들려다 원자가 위반도 났다. 1급→이민, 2급→엔아민으로 분리했다."},
        {"d": _EXC, "k": "fix", "t": "산물 폭발 억제",
         "b": "암로디핀에서 1,113종이 나왔고 657종(59%)이 에스터교환이었다 "
              "(에스터 2개 × 폴리올 OH 다수). 규칙·쌍마다 최대 3개로 위치이성질체를 접는다."},
    ],
    "p19_excipient_library_ext": [
        {"d": _EXC, "k": "new", "t": "신규 — 부형제 76종 → 143종",
         "b": "구조 유효 46→143종, 반응 조각 보유 39→102종, "
              "<b>예측에 참여 못 하던 항목 25종 → 0종</b>. "
              "원본 <code>p70</code>은 건드리지 않고 <code>apply_extensions()</code>로 덮어쓴다."},
        {"d": _EXC, "k": "fix", "t": "죽어 있던 25종을 살림",
         "b": "메타중아황산나트륨이 여기 있었다. 아황산염은 알데하이드·케톤 원료와 "
              "부가물을 만드는 <b>문헌에 잘 정리된 배합금기</b>인데 «이상 없음»으로 "
              "보고되고 있었다. 젤라틴 · 레시틴 · 소르빈산칼륨도 마찬가지."},
        {"d": "2026-08-26", "k": "fix", "t": "PEG 4000 과산화물 조각 누락",
         "b": "아미오다론 시험 중 발견. 같은 계열 부형제는 모두 <code>OO</code>를 "
              "가지고 있는데 PEG 4000만 빠져 있었다. <code>PROBE_PATCH</code>로 보강."},
        {"d": "2026-08-26", "k": "fix", "t": "구조 결함 5건",
         "b": "스테아르산마그네슘 SMILES가 사슬 전체를 대괄호로 묶은 문법 오류였다 "
              "(<code>[CCCCCCCCCCCCCCCCCC(=O)[O-]]</code>). 빈 구조 4건도 함께 수정."},
        {"d": "2026-08-26", "k": "note", "t": "고분자의 smiles는 대표 조각",
         "b": "실제 고분자 구조가 아니다. 몰비 계산에 쓰지 않고 화면에 «정의 없음»으로 표시한다."},
    ],

    # ---------------------------------------------------- 3x
    "p30_reaction_engine": [
        {"d": _P30, "k": "fix", "t": "3급 아민 니트로소화 반응식 재작성",
         "b": "생성물의 <code>N=O</code>가 <b>질소가 아니라 탄소에</b> 붙고 있었다. "
              "N-니트로소아민이 아니라 C-니트로소 화합물이 나와 하류 필터 "
              "<code>[#7]-[#7]=[#8]</code>가 전부 걸러냈다 — 3급 아민 원료의 NDSRI가 "
              "항상 0건이었던 원인. <b>니트로소 탈알킬화</b>로 다시 썼다 "
              "(R₂N–CH₂R′ + HNO₂ → R₂N–N=O + R′CHO). "
              "아미오다론 → NDSRI + NDEA, 디펜히드라민 → NDSRI + NDMA."},
        {"d": _P30, "k": "fix", "t": "아민 탐지 SMARTS 방향족 오탐 제거",
         "b": "<code>[#7]</code>은 방향족 질소(소문자 n)도 잡는다. 카페인이 3급 아민 3개로, "
              "피롤·인돌이 2급 아민으로 세어지고 있었다. <code>[NX3]</code>로 한정해 "
              "<b>참양성 손실 없이 위양성 5건</b>을 제거했다. 아테놀올·프로프라놀올 결과는 그대로."},
        {"d": _P30, "k": "fix", "t": "reagent_ndsri_user 규칙 폐기",
         "b": "결함이 셋이었다. ① 전하가 틀림 — 아질산 이온의 <code>[O-]</code>만 떼고 "
              "질소의 +1을 물려받아 상대 음이온 없는 <code>[N+]=O</code>가 생성됨 "
              "② <code>[H:3]</code> 매핑 탓에 명시적 수소가 남아 DB 표준 SMILES와 대조되지 않음 "
              "③ 아마이드 배제 조건이 없어 아테놀올의 1차 아마이드까지 니트로소화 "
              "(<code>C(=O)N([H])[N+]=O</code>). "
              "하류 필터 <code>[#7]-[#7]=[#8]</code>는 전하를 보지 않아 이 구조들이 "
              "<b>화면까지 올라오고 있었다</b> (아테놀올+PEG4000 11건 중 3건). "
              "같은 산물의 올바른 형태는 <code>reagent_ndsri_secondary</code>가 이미 만들므로 "
              "제거했다. 원료 10종 × 부형제 8종 전수 대조에서 "
              "<b>정상 산물 손실 0건</b>, 총 산물 171 → 163건, 위험도 밴드 변화 없음. "
              "SMIRKS 문자열은 <code>RETIRED_…</code>로 기록만 남겼다."},
    ],

    # ---------------------------------------------------- 6x
    "p63_deepchem_predictor": [
        {"d": "2026-08-12", "k": "fix", "t": "deepchem · tensorflow 의존 제거",
         "b": "설치가 까다롭고 버전 충돌이 잦아 RDKit ECFP(r=2 · 2048bit) + "
              "scikit-learn RandomForest로 바꿨다. 파일명은 호출부 호환을 위해 유지. "
              "원본은 <code>deepchem_predictor.py.backup</code>."},
        {"d": "2026-08-12", "k": "note", "t": "Ames 데이터로 재학습",
         "b": "외부 검증 n=1,440 (양성 785 · 음성 655)에서 민감도 86.0% · 특이도 79.4%."},
    ],
    "p62_toxicity_alerts": [
        {"d": "2026-08-21", "k": "note", "t": "특이도 상한의 원인",
         "b": "Ashby-Tennant 규칙 기반. 외부 검증에서 민감도 70.8% · 특이도 52.7%이며 "
              "위양성 310건이 앙상블 특이도 상한을 만든다."},
    ],
    "p66_ich_m7": [
        {"d": "2026-08-21", "k": "note", "t": "상향 전용 합의",
         "b": "세 층 중 하나라도 위험하다고 하면 올린다. 내리지 않는다. "
              "통계적으로 상향된 161종 중 146종이 실제 Ames 양성이었다(정밀도 90.7%). "
              "하향을 허용하는 변형은 <code>engines/p67</code>."},
    ],

    # ---------------------------------------------------- 7x
    "p70_excipient_master_library": [
        {"d": _EXC, "k": "note", "t": "이 파일은 원본 그대로",
         "b": "확장은 <code>p19_excipient_library_ext.apply_extensions()</code>가 "
              "런타임에 덮어쓴다. 원본 값이 유효하면 덮어쓰지 않는다."},
    ],
    "p72_excipient_interaction_engine": [
        {"d": _EXC, "k": "fix", "t": "반응 규칙을 p18에서 불러오도록 변경",
         "b": "내장 2종 대신 <code>p18.build_reaction_rules()</code> 14종을 쓴다. "
              "실패하면 기존 2종(<code>_FALLBACK_RULES</code>)으로 떨어진다. "
              "라이브러리도 <code>p19</code> 확장을 적용한다. 백업 <code>_backup_excipient_ext/</code>."},
        {"d": "2026-08-26", "k": "new", "t": "⭐ p73을 흡수 — 니트로소화 경로 단일화",
         "b": "니트로소화 엔진이 두 벌이었고 서로 다른 기준으로 판단해 "
              "<b>«한쪽이 놓치는 것을 다른 쪽이 잡는»</b> 상태였다. "
              "p73은 API 게이트(<code>[NX3;H1,H2]</code>)가 3급 아민을 막고, "
              "이쪽은 부형제 게이트(한글 경고문 키워드)가 크로스포비돈·MCC를 막았다. "
              "둘이 <code>seen_smiles</code>를 공유해 먼저 도는 쪽이 산물을 가져가기까지 했다. "
              "이제 <code>_evaluate_nitrosamine_risk</code> 하나뿐이다. "
              "아질산염 매트릭스 7종도 이 파일로 옮겼다. 542줄 → 694줄."},
        {"d": "2026-08-26", "k": "fix", "t": "부형제 게이트를 프로파일 우선으로",
         "b": "한글 경고문 키워드 대신 아질산염 잔류 프로파일을 먼저 본다. "
              "<code>has_nitrite_risk=False</code>는 «자료 없음»이 아니라 "
              "«평가했고 무시할 수준»이므로 차단한다(유당·스테아르산마그네슘). "
              "프로파일이 없으면 반응 조각 → 키워드 순으로 폴백하되 "
              "«⚠ 실측 프로파일 없음 — trace 가정»으로 화면에 명시한다. "
              "<b>크로스포비돈(0.1–5.0 ppm)과 MCC가 처음으로 니트로소화 경로에 들어왔다.</b>"},
        {"d": "2026-08-26", "k": "note", "t": "병합 결과",
         "b": "원료 8종 × 부형제 143종 전수 대조에서 <b>1,432 → 1,451건 (+19), 손실 0건</b>. "
              "3급 아민 원료(아미오다론·디펜히드라민·베라파밀)에서 니트로소화가 도는 부형제가 "
              "7종 → 9종으로 늘었다. 위험도 순위도 잔류 ppm을 따르게 됐다 — "
              "크로스포비돈 88 &gt; PEG 4000 85 &gt; MCC 79."},
    ],

    # ---------------------------------------------------- 8x
    "p81_report_export": [
        {"d": _FLASK, "k": "fix", "t": "한글 PDF가 Windows에서만 되던 문제",
         "b": "<code>_windows_fonts()</code>가 <code>C:\\Windows\\Fonts</code>만 찾았다. "
              "폰트를 못 찾으면 fpdf가 Helvetica로 떨어지는데 Helvetica는 latin-1만 지원하므로 "
              "<b>한글이 한 글자라도 있으면 PDF 생성이 실패</b>한다. "
              "리눅스(나눔고딕 · Noto Sans CJK)와 맥 경로를 추가했다. Windows 동작은 그대로."},
    ],
    "p88_excipient_report_export": [
        {"d": _EXC, "k": "open", "t": "본문이 아직 부분적으로 한글",
         "b": "화면은 EN/KO가 되지만 이 보고서 본문은 영문화가 끝나지 않았다. "
              "<code>p86_degradation_report_export</code>도 같은 상태."},
    ],
    "p89_drug_db_report_export": [
        {"d": "2026-08-26", "k": "fix", "t": "Flask 판에 다시 연결",
         "b": "Flask 전환 때 <b>이 모듈을 부르는 곳이 한 군데도 없어졌다.</b> "
              "Streamlit 판에서는 <code>p85_report_download_ui</code>가 이어 주고 있었는데 "
              "그 모듈이 화면과 함께 <code>_archive</code>로 가면서 DB 검색 리포트가 통째로 끊겼다. "
              "이 지도를 다시 만들면서 «아무도 부르지 않는 모듈»로 드러났다. "
              "<code>reports.drugdb_preview</code> · <code>reports.drugdb_download</code> "
              "라우트를 붙이고 검색 화면과 보고서 보관함에 버튼을 넣었다. "
              "TXT · Word · PDF · Excel · JSON 5종 모두 생성 확인."},
        {"d": "2026-08-26", "k": "fix", "t": "streamlit 의존 제거",
         "b": "<code>render_drug_db_report_downloads()</code>(L384~481)가 함수 안에서 "
              "<code>import streamlit as st</code>를 하고 있었다. Flask 라우트가 같은 일을 하므로 "
              "통째로 걷어냈다(481줄 → 382줄). "
              "<b>이제 코드베이스에 streamlit을 참조하는 모듈이 없다.</b>"},
        {"d": "2026-08-26", "k": "note", "t": "화면 표와 열 이름이 다르다",
         "b": "화면은 소문자 키(<code>name · smiles · m7_class</code>)를 쓰고 p89는 "
              "표시용 열 이름(<code>Name · SMILES · ICH M7</code>)을 기대한다. "
              "<code>webapp/services/search.drug_db_report_context()</code>가 사이에서 옮긴다. "
              "<code>공정 안정성</code>은 p50에서 가져오는데 값에 그림문자(✅)가 섞여 있어 "
              "PDF 폰트에 글리프가 없으므로 걷어낸 뒤 넣는다."},
    ],
    "p86_degradation_report_export": [
        {"d": _EXC, "k": "open", "t": "본문이 아직 부분적으로 한글",
         "b": "영문 보고서 요청 시 한글이 섞여 나온다."},
    ],

    # ---------------------------------------------------- 스크립트
    "s7_run_case_study": [
        {"d": "2026-08-26", "k": "fix", "t": "⭐ 라니티딘·니자티딘에서 NDMA가 나온다",
         "b": "p30 3급 아민 규칙을 고친 뒤 재실행했더니 "
              "<code>CN(C)N=O</code>(<b>NDMA</b> · Class 1 · p=0.96)가 새로 잡혔다. "
              "<b>이전 판은 NDMA를 놓치고 있었다.</b> 전구체는 다이메틸아미노기 — 즉 3급 아민인데 "
              "규칙이 C-니트로소를 만들어 하류 필터가 전부 걸러냈기 때문이다 "
              "(구 로그: «제외 — 전하 표기 오류 2 · N-니트로사민 아님 2»). "
              "보고된 기전(3급 아민 → 불안정 니트로소암모늄 → 탈알킬화 → DMA → NDMA)과 일치한다."},
        {"d": "2026-08-26", "k": "note", "t": "니트로엔아민 산물 2종이 빠졌다",
         "b": "이전 산물은 라니티딘의 <b>니트로엔아민 질소</b>에 니트로소기를 붙인 것이었다. "
              "나이트로기와 공액된 push-pull 계라 아마이드에 가까울 만큼 친핵성이 약하다. "
              "새 <code>SECONDARY_AMINE_SMARTS</code>가 치환기를 <code>[CX4,c]</code>로 한정하면서 "
              "sp2 탄소에 붙은 질소가 제외됐다."},
        {"d": "2026-08-26", "k": "open", "t": "논문 본문 수정 필요",
         "b": "사례연구 표에 Class 1이 2종 생긴다. 구버전은 "
              "<code>논문자료/03_사례연구/_구버전_20260821/</code> 에 보관."},
    ],
    "s8_run_excipient_case_study": [
        {"d": "2026-08-26", "k": "note", "t": "재실행 — 15건 → 44건",
         "b": "8/25 부형제 확장·반응규칙 14종까지 누적된 결과다. 사라진 8건은 전부 원인이 확인된다 — "
              "아미드 N 이민 4건(Maillard 규칙 수정) · 니트로엔아민 N 니트로소화 4건(2급 SMARTS 강화). "
              "설명되지 않는 손실은 없다."},
    ],
    "s9_run_binary_compatibility_study": [
        {"d": "2026-08-26", "k": "note", "t": "재실행 — 가려진 쌍 0쌍",
         "b": "5종 패널에서는 이제 가려지는 조합이 없다. 143종 전수로 재면 "
              "634 → 1,021건 (+387), 원료마다 48~49종이 가려진다. "
              "표본이 작으면 충돌 여부가 우연에 좌우된다는 뜻이므로 "
              "논문에는 전수 수치를 쓰는 편이 낫다."},
        {"d": _EXC, "k": "new", "t": "신규 — 이원 배합 사례연구",
         "b": "<code>s8</code>의 사례를 1:1 이원 모드로 다시 돌리고 "
              "<code>binary_vs_matrix_comparison.csv</code>를 만든다. 논문 표 5 갱신용."},
    ],
    "s10_build_module_map": [
        {"d": "2026-08-26", "k": "new", "t": "신규 — 이 지도를 만드는 스크립트",
         "b": "코드베이스를 AST로 훑어 HTML 안의 DATA 블록만 교체한다. "
              "역할 한 줄과 변경 이력은 이 파일의 <code>ROLES</code> · <code>HIST</code>에 적는다."},
    ],

    # ---------------------------------------------------- webapp
    "webapp/__init__": [
        {"d": _FLASK, "k": "new", "t": "신규 — Flask 앱 팩토리",
         "b": "블루프린트 등록과 설정을 담당한다. Streamlit 판에는 대응물이 없다."},
    ],
    "webapp/state": [
        {"d": _FLASK, "k": "new", "t": "신규 — st.session_state 대체",
         "b": "서버 메모리에 세션을 둔다. TTL 8시간. "
              "<code>bucket() · get/put/pop · display_prefs() · canvas()</code>. "
              "<b>프로그램을 껐다 켜면 캔버스와 결과가 사라진다.</b>"},
    ],
    "webapp/i18n": [
        {"d": _I18N, "k": "new", "t": "신규 — 화면 번역 사전 461키",
         "b": "영어 기본, 오른쪽 위 <code>EN | KO</code> 버튼으로 전환. "
              "번역이 없는 키는 키 이름을 그대로 보여주고 로그에 경고를 남긴다 — "
              "누락을 조용히 넘기지 않기 위해서."},
        {"d": _I18N, "k": "note", "t": "엔진 데이터는 건드리지 않았다",
         "b": "p42 · p43 · p68 · p72 등에는 한글 라벨이 박혀 있고 값 자체가 식별자로 쓰이는 "
              "곳이 있어 바꾸면 로직이 깨진다. 화면에 보일 때만 <code>opt()</code>로 변환한다."},
        {"d": _P30, "k": "fix", "t": "p30 새 규칙 라벨 추가",
         "b": "«니트로소 탈알킬화» 대응 영문을 넣었다. 구 라벨도 남겨 뒀다."},
    ],
    "webapp/excipient_terms": [
        {"d": _I18N, "k": "new", "t": "신규 — 부형제 영문명 · 주의문구",
         "b": "처음 76종으로 시작해 라이브러리 확장에 맞춰 143종으로 늘렸다. "
              "<code>COLUMN_EN · NITRITE_NOTE_EN · name_en() · coverage()</code>."},
    ],
    "webapp/services/chem": [
        {"d": _FLASK, "k": "fix", "t": "예외를 삼키던 판정 함수",
         "b": "<code>classify()</code>가 예외를 삼키고 <code>None</code>만 돌려줘 "
              "화면에 «판정 실패»만 뜨고 이유를 알 수 없었다. "
              "<code>(result, error)</code> 튜플로 바꾸자마자 시험 중 DB 접근 실패가 드러났다."},
    ],
    "webapp/services/search": [
        {"d": "2026-08-26", "k": "new", "t": "DB 검색 리포트 연결",
         "b": "<code>drug_db_report_context()</code> · <code>drug_db_report_text()</code>를 추가해 "
              "<code>p89</code>를 되살렸다. 화면 표(소문자 키)를 p89의 표시용 열 이름으로 옮기고, "
              "파일명에서 의약품 이름을 뽑아내며(<code>Metformin_related.csv</code> → "
              "<code>Metformin</code>), PDF 폰트에 없는 그림문자를 걷어낸다."},
        {"d": _FLASK, "k": "new", "t": "신규 — p11 검색 이력을 대체",
         "b": "<code>p11_structure_search_history</code>는 <code>st.session_state</code>를 "
              "직접 조작하므로 옮길 수 없었다. 같은 기능을 새로 구현했다."},
        {"d": _FLASK, "k": "fix", "t": "반환 타입 불일치 정리",
         "b": "<code>substructure_search · similarity_search_top_n</code>은 튜플을, "
              "<code>load_related_table</code>은 DataFrame을 돌려준다. "
              "<code>resolve_data_path(stored, file_name)</code>은 인자가 2개다."},
    ],
    "webapp/services/excipient": [
        {"d": "2026-08-26", "k": "fix", "t": "카테고리 필터가 아예 작동하지 않았음",
         "b": "<code>filter_excipient_keys</code>는 키워드 전용 인자를 받는데 위치 인자로 "
              "부르고 있었다. <code>TypeError</code>가 나고 그걸 "
              "<code>except Exception: pass</code>가 삼켰다. 서비스에서 다시 구현."},
        {"d": "2026-08-26", "k": "fix", "t": "143종 중 일부만 보이던 문제",
         "b": "420px 스크롤 상자에 묻혀 있었다. 카테고리별 아코디언 + 필터 칩 + "
              "펼치기/접기/전체선택으로 다시 만들었다. 카탈로그 수가 76에 고정돼 있던 것도 "
              "(<code>opensource_catalog_stats()</code>에서 오던 값) 실제 선택지로 다시 센다."},
    ],
    "webapp/services/process": [
        {"d": _FLASK, "k": "fix", "t": "공정 xlsx가 0바이트로 받아지던 문제",
         "b": "경로 B에 기록이 없으면 빈 파일이 생성됐다. 가드를 넣고 버튼을 비활성 처리했다."},
    ],
}

# 이력이 없는 모듈은 지도에서 «변경 기록 없음»으로 표시된다.

# ---------------------------------------------------------------------------
# 수집
# ---------------------------------------------------------------------------

SKIP_DIRS = {
    "_archive", "_backup_before_rename", "_backup_excipient_ext",
    "RxnMetrics", "__pycache__", ".git", "venv", "node_modules",
    "data", "cache", "models", "논문자료", "md파일", "templates",
}

STDLIB = set(getattr(sys, "stdlib_module_names", ()))
STDLIB |= {"typing", "dataclasses", "enum", "pathlib", "collections"}


def _iter_module_files() -> list[tuple[str, str]]:
    """(module_id, abs_path) 목록."""
    found: list[tuple[str, str]] = []

    for name in sorted(os.listdir(ROOT)):
        if name.endswith(".py") and re.match(r"^[ps]\d", name):
            found.append((name[:-3], os.path.join(ROOT, name)))

    eng = os.path.join(ROOT, "engines")
    if os.path.isdir(eng):
        for name in sorted(os.listdir(eng)):
            if name.endswith(".py") and re.match(r"^p\d", name):
                found.append((name[:-3], os.path.join(eng, name)))

    web = os.path.join(ROOT, "webapp")
    for dirpath, dirnames, filenames in os.walk(web):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for name in sorted(filenames):
            if not name.endswith(".py"):
                continue
            full = os.path.join(dirpath, name)
            rel = os.path.relpath(full, ROOT).replace(os.sep, "/")[:-3]
            if os.path.getsize(full) < 60 and name == "__init__.py":
                continue                                   # 빈 패키지 마커는 뺀다
            found.append((rel, full))

    return found


def _analyse(path: str) -> dict:
    src = open(path, encoding="utf-8").read()
    tree = ast.parse(src)

    pub: list[str] = []
    const: list[str] = []
    imports: list[tuple[str, list[str]]] = []

    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if not node.name.startswith("_"):
                pub.append(node.name)
        elif isinstance(node, ast.Assign):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and re.match(r"^[A-Z][A-Z0-9_]*$", tgt.id):
                    const.append(tgt.id)
        elif isinstance(node, ast.AnnAssign):
            if isinstance(node.target, ast.Name) and re.match(r"^[A-Z][A-Z0-9_]*$", node.target.id):
                const.append(node.target.id)

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                imports.append((a.name, [a.asname or a.name.split(".")[0]]))
        elif isinstance(node, ast.ImportFrom):
            if node.level:                                  # 상대 import
                imports.append(("." * node.level + (node.module or ""),
                                [a.name for a in node.names]))
            elif node.module:
                imports.append((node.module, [a.name for a in node.names]))

    return {
        "src": src,
        "loc": len(src.splitlines()),
        "bytes": len(src.encode("utf-8")),
        "pub": pub,
        "const": const,
        "imports": imports,
    }


def _band_of(mid: str) -> str:
    if mid.startswith("webapp/views"):
        return "bwv"
    if mid.startswith("webapp/"):
        return "bws"
    if mid.startswith("s"):
        return "bs"
    m = re.match(r"^p(\d)", mid)
    return "b" + m.group(1) if m else "b9"


def _short_of(mid: str) -> str:
    if mid.startswith("webapp/views/"):
        return mid.split("/")[-1]
    if mid.startswith("webapp/services/"):
        return mid.split("/")[-1]
    if mid.startswith("webapp/"):
        return mid.split("/")[-1] if mid != "webapp/__init__" else "app factory"
    m = re.match(r"^([ps]\d+)_", mid)
    return m.group(1) if m else mid


def _pkg_of(mid: str, path: str) -> str:
    rel = os.path.relpath(path, ROOT).replace(os.sep, "/")
    return "/".join(rel.split("/")[:-1]) or "root"


def _resolve(target: str, mid: str, ids: set[str]) -> str | None:
    """import 대상 문자열을 모듈 id로."""
    if target.startswith("."):                              # 상대 import
        base = mid.rsplit("/", 1)[0]
        stripped = target.lstrip(".")
        up = len(target) - len(stripped) - 1
        parts = base.split("/")
        if up:
            parts = parts[:-up] if up <= len(parts) else []
        cand = "/".join([p for p in parts if p] + ([stripped.replace(".", "/")] if stripped else []))
        for c in (cand, cand + "/__init__"):
            if c in ids:
                return c
        return None
    if target in ids:
        return target
    dotted = target.replace(".", "/")
    # engines/ 안의 모듈은 id에 폴더명을 붙이지 않았다 (구버전 지도와 맞추기 위해).
    # 실제 코드는 `from engines.p33_... import ...` 형태로 부르므로 접두사를 떼어 준다.
    cands = [dotted, dotted + "/__init__"]
    if dotted.startswith("engines/"):
        bare = dotted[len("engines/"):]
        cands += [bare, bare + "/__init__"]
    for c in cands:
        if c in ids:
            return c
    return None


def build() -> dict:
    files = _iter_module_files()
    ids = {mid for mid, _ in files}
    info = {mid: _analyse(path) for mid, path in files}
    paths = dict(files)

    nodes = []
    missing_role = []
    for mid, path in files:
        a = info[mid]
        ext = sorted({
            t.split(".")[0] for t, _ in a["imports"]
            if not t.startswith(".")
            and _resolve(t, mid, ids) is None
            and t.split(".")[0] not in STDLIB
        })
        role = ROLES.get(mid)
        if role is None:
            role = "(설명 없음)"
            missing_role.append(mid)
        doc = ast.get_docstring(ast.parse(a["src"])) or ""
        nodes.append({
            "id": mid,
            "pkg": _pkg_of(mid, path),
            "bytes": a["bytes"],
            "loc": a["loc"],
            "pub": a["pub"],
            "const": a["const"],
            "ext": ext,
            "role": role,
            "band": _band_of(mid),
            "short": _short_of(mid),
            "name": mid.split("/")[-1],
            "desc": " ".join(doc.split())[:400],
            "hist": HIST.get(mid, []),
        })

    edge_map: dict[tuple[str, str], set[str]] = {}
    for mid, _ in files:
        for target, syms in info[mid]["imports"]:
            pkg_dest = _resolve(target, mid, ids)
            for sym in syms:
                # `from ..services import single` 처럼 심볼 자체가 모듈인 경우가 있다.
                # 패키지 __init__ 이 아니라 그 모듈로 선을 그어야 한다.
                sep = "" if target.endswith(".") else "."
                sub = _resolve(target + sep + sym, mid, ids)
                dest = sub or pkg_dest
                if dest and dest != mid:
                    edge_map.setdefault((mid, dest), set()).add(sym)

    edges = [
        {"s": s, "t": t, "syms": sorted(v)[:20], "n": len(v)}
        for (s, t), v in sorted(edge_map.items())
    ]

    if missing_role:
        print("⚠ ROLES 에 없는 모듈 — 설명을 적어 주세요:")
        for m in missing_role:
            print("   ", m)

    stale = [k for k in list(ROLES) + list(HIST) if k not in ids and not k.endswith(".py")]
    if stale:
        print("⚠ 코드에 없는데 ROLES/HIST 에 남아 있는 항목:")
        for m in stale:
            print("   ", m)

    used = {n["band"] for n in nodes}
    bands = [b for b in BANDS if b["id"] in used]

    return {"nodes": nodes, "edges": edges, "bands": bands}


def inject(data: dict) -> None:
    html = open(HTML, encoding="utf-8").read()
    key = "const DATA = "
    i = html.index(key) + len(key)
    depth = 0
    for k in range(i, len(html)):
        if html[k] == "{":
            depth += 1
        elif html[k] == "}":
            depth -= 1
            if depth == 0:
                j = k + 1
                break
    payload = json.dumps(data, ensure_ascii=False, separators=(", ", ": "))
    open(HTML, "w", encoding="utf-8").write(html[:i] + payload + html[j:])


if __name__ == "__main__":
    d = build()
    inject(d)
    nb = {}
    for n in d["nodes"]:
        nb[n["band"]] = nb.get(n["band"], 0) + 1
    print(f"\n노드 {len(d['nodes'])}개 · 엣지 {len(d['edges'])}개 · 밴드 {len(d['bands'])}개")
    for b in d["bands"]:
        print(f"   [{b['id']:4s}] {nb.get(b['id'],0):3d}개  {b['label']}")
    h = sum(1 for n in d["nodes"] if n["hist"])
    print(f"변경 이력이 붙은 모듈 {h}개 / 전체 {len(d['nodes'])}개")
    print(f"→ {HTML}")
