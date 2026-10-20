"""
Ashby-Tennant 유전독성 구조 알러트 탐지 (SMARTS) 및 하이라이트 좌표 생성.

학계·ICH M7 Option 4 (Expert Rule-based)에서 인용되는 대표 알러트 패턴.

THIRD-PARTY ATTRIBUTION — Toxtree (GPL-2.0)
-------------------------------------------
Some patterns in this module derive from the Benigni-Bossa rule base as
implemented in Toxtree 3.1.0, (c) Ideaconsult Ltd. and contributors,
licensed under GPL-2.0. Each is marked at its point of definition with one
of three tags:

    Attribution: verbatim from Toxtree 3.1.0 <file>
        the SMARTS string was transcribed unchanged
    Attribution: adapted from Toxtree 3.1.0 <file>
        the Toxtree pattern was narrowed or rewritten for this engine
    Attribution: independent reimplementation of Toxtree 3.1.0 <file>
        only the decision logic was reimplemented; no code or pattern copied

Grep for "Attribution:" to list every one of them. No part of the Toxtree
distribution is included in or redistributed with this project, and Toxtree
was run as a benchmark only, which is use rather than distribution. See the
"Third-party material and attribution" section of README.md.

The remaining patterns are written from the primary literature (Ashby &
Tennant; Benigni & Bossa, Chem. Rev. 2011) or from ICH M7(R2) Section 7.5,
and carry no Toxtree attribution.
"""

from __future__ import annotations

from rdkit import Chem

# 여러 알러트 목록에서 공유하는 SMARTS — 문자열을 그대로 재사용해야
# SMARTS 기준 중복 제거가 동작합니다 (같은 화학적 실체가 두 번 계상되지 않도록).
# 2026-09-03 수정: 아민 질소 자신을 배제하던 재귀 조건 `[#7;!$(N[*]=O)]` 때문에
# NDMA(CN(C)N=O)가 매치되지 않던 것을 수정 (니트로소 질소만 배제).
NITROSAMINE_SMARTS = "[#7;!$([#7]=O);!$([#7+][O-])]-[#7X2]=[#8]"
# 2026-09-03 수정: `[N;$(N-c1ccccc1)]`은 페닐에 붙은 모든 질소(아마이드·니트로)를
# 잡았습니다. 1·2급 아민만 남기고 아마이드/설폰아마이드/히드라진/N-O를 배제.
AROMATIC_AMINE_SMARTS = "[NX3;H2,H1;!$(N[C,S,P]=[O,S]);!$(N-[N,O])]-c"
# 2026-09-03 수정: 방향족 니트로와 이중 계상되지 않도록 방향족에 붙은 니트로는 제외.
NITRO_AROMATIC_SMARTS = "[c,n][NX3+](=O)[O-]"
NITRO_ALIPHATIC_SMARTS = "[$([NX3+](=O)[O-]);!$([NX3+](=O)([O-])[c,n])]"
# 2026-09-14 수정: `C=CC(=O)` 는 너무 헐거웠다. 학습셋에서 431건을 잡으면서
# PPV 48.0% — 전체 양성 비율 54.5% 보다 낮아 정보가 없는 경보였고, 위양성
# 224건의 최대 단일 원인이었다.
#
# Benigni/Bossa SA10_gen 의 패턴으로 교체한다 (Toxtree 3.1.0 소스,
# mutant/rules/SA10_gen.java 의 AB_UNSATURATED_CARBONYLS_SMARTS 문자열 그대로).
# Attribution: verbatim from Toxtree 3.1.0 mutant/rules/SA10_gen.java
#   (AB_UNSATURATED_CARBONYLS_SMARTS). GPL-2.0, (c) Ideaconsult Ltd.
#   This alert is currently disabled — see the note at its entry below.
# 배제되는 것 —
#   · 퀴논          `!$(C1(=O)C=CC(=O)C=C1)` — 별도 경보가 이미 담당
#   · 카복실산·에스터  카보닐의 다른 치환기가 산소면 제외 `[!O;...]`
#   · 방향족 고리의 이중결합  `!:;=` 로 방향족 결합을 제외
#
# ⚠️ SA10 은 «β 탄소의 치환기에 탄소가 6개 이상이면 제외»하는 규칙을 Java 쪽에서
#    따로 구현한다. SMARTS 로는 표현할 수 없어 여기에는 반영되지 않았다.
#    그만큼 Toxtree 보다 넓게 잡힌다.
MICHAEL_ACCEPTOR_SMARTS = (
    "[!a,#1;!$(C1(=O)C=CC(=O)C=C1)][#6]([!a,#1;!$(C1(=O)C=CC(=O)C=C1)])"
    "!:;=[#6][#6](=O)[!O;!$([#6]1:,=[#6][#6](=O)[#6]:,=[#6][#6](=O)1)]"
)

# Ashby-Tennant 핵심 알러트 — SMARTS (규제·문헌 참고)
ASHBY_TENNANT_ALERTS: tuple[tuple[str, str], ...] = (
    # 사용자 요청 핵심 패턴
    # 2026-09-03 수정: `[C,c][Cl,Br,I]`의 소문자 c 때문에 클로르프로마진·로사르탄
    # 등 비반응성 아릴 할라이드가 전부 POSITIVE로 잡혔습니다. sp3 탄소로 한정하고
    # 활성화된 아릴 할라이드는 아래 별도 항목으로 분리.
    ("Alkylating agent (알킬화제, sp3 halide)", "[CX4;!$(C=[O,N,S])][Cl,Br,I]"),
    (
        "Activated aryl halide (활성화 아릴 할라이드)",
        "[F,Cl,Br,I]-[c;$(c1ccccc1[N+](=O)[O-]),$(c1ccc([N+](=O)[O-])cc1),"
        "$(c1ncccn1),$(c1ccncc1),$(c1ncccc1)]",
    ),
    ("Aromatic amine (방향족 아민)", AROMATIC_AMINE_SMARTS),
    ("Nitro group (aromatic)", NITRO_AROMATIC_SMARTS),
    ("Nitro group (aliphatic)", NITRO_ALIPHATIC_SMARTS),
    # 확장 Ashby-Tennant 대표 알러트
    # 2026-09-03 수정: 이름만 다르고 화학적으로 동일하던 "N-nitroso (N-N=O)" 항목 삭제
    # (하나의 니트로소기가 두 건으로 계상되던 원인).
    ("N-nitroso (니트로사민)", NITROSAMINE_SMARTS),
    ("Aromatic N-oxide", "[n+][O-]"),
    # 2026-09-14 수정: 에폭사이드만 보던 `[OX2r3]`을 3원환 헤테로고리로 확장.
    # Benigni/Bossa SA7_gen 은 에폭사이드와 아지리딘을 하나의 패턴 `C1[O,N]C1`로
    # 처리한다. 학습셋 결손 분석에서 아지리딘 계열 40건이 누락되고 있었다.
    # 원 패턴은 C 에 X4·H 제약이 없어 매우 넓게 매칭되므로, 고리 원자로 한정하고
    # 질소는 아마이드성 질소를 배제해 위양성을 줄였다.
    # Attribution: adapted from Toxtree 3.1.0 mutant/rules/SA7_gen.java
    #   (one C1[O,N]C1 pattern; narrowed here to ring atoms, amide N excluded).
    ("Epoxide / aziridine (에폭사이드·아지리딘)", "[CX4;r3][OX2,NX3;r3;!$([NX3]C=[O,S])][CX4;r3]"),
    # 2026-09-03 수정: "Alkyl halide [Cl,Br,I][CH2;!$(C=O)]"는 위 sp3 알킬화제 패턴의
    # 진부분집합이라 같은 결합을 두 번 계상하므로 삭제.
    ("Aldehyde (알데하이드)", "[CH1](=O)[#6]"),
    # 2026-09-14 수정: `O=C1C=CC(=O)C=C1` 은 이중결합만 받아 방향족으로 표기된
    # 고리를 놓쳤다. 학습셋에서 15건밖에 잡지 못했고 안트라퀴논 계열 59건이
    # 결손으로 남았다. Benigni/Bossa SA12_gen 의 두 패턴으로 교체한다.
    # `=,:` 가 이중결합과 방향족을 모두 받으므로 융합 고리의 퀴논도 매칭된다.
    # Attribution: adapted from Toxtree 3.1.0 mutant/rules/SA12_gen.java
    #   (para and ortho quinone patterns; bond primitives relaxed to `=,:` here
    #   so that quinones drawn as aromatic rings also match).
    ("Quinone (퀴논, para)", "O=[#6]1[#6]=,:[#6][#6](=O)[#6]=,:[#6]1"),
    ("Quinone (퀴논, ortho)", "O=[#6]1[#6]=,:[#6][#6]=,:[#6][#6]1(=O)"),
    ("Hydrazine (히드라진)", "[NX3][NX3]"),
    # 2026-09-03 수정: `[CH2OH]`는 대괄호 안이 AND 결합이라 «탄소이면서 산소»인
    # 원자를 요구해 아무것도 매치하지 못했습니다 (알러트가 상시 무음).
    ("Methylol (N-methylol)", "[OX2H]-[CH2]-[NX3]"),
    ("Azo aromatic (아조 방향족)", "[c][N;X2]=[N;X2][c]"),
    # 2026-09-14 제외: Michael acceptor.
    #
    # 학습셋 근거 —
    #   `C=CC(=O)` (헐거운 원 패턴)   431건 · PPV 48.0%
    #   SA10 패턴으로 좁힌 뒤          60건 · PPV 35.0%
    # 양성 비율 54.5% 를 두 경우 모두 밑돈다. 좁힐수록 PPV 가 오히려 떨어졌는데,
    # 이는 맞는 쪽부터 먼저 탈락했다는 뜻이다.
    #
    # 원인 추정 — SA10_gen 은 «β 탄소의 치환기에 탄소가 6개 이상이면 제외»를
    # Java 쪽에서 따로 구현한다. SMARTS 로 옮길 수 없어 그 거름망이 빠졌고,
    # 결과적으로 Toxtree 수준의 정밀도에 도달하지 못했다.
    #
    # ⚠️ 마이클 수용체는 고전적 구조경보이므로 제외 사실을 논문에 명시할 것.
    #    β 탄소 규칙을 파이썬으로 구현하면 되살릴 수 있다.
    # ("Michael acceptor (α,β-불포화 카르보닐)", MICHAEL_ACCEPTOR_SMARTS),
    # 2026-09-03 수정: ICH M7 7.5절 cohort of concern(aflatoxin-like·N-nitroso·
    # alkyl-azoxy) 중 aflatoxin·azoxy 알러트가 없어 아플라톡신 B1이 «알러트 0건»
    # → Class 5(비변이원성)로 떨어지던 문제. 모체 구조에서 잡히는 알러트를 추가.
    (
        "Aflatoxin-like (coumarin-fused bis-furan)",
        "[CX4;R]1[#6;R]~[#6;R][OX2;R][CX4H1;R]1[OX2;R][c;R]",
    ),
    # Attribution: related to Toxtree 3.1.0 mutant/rules/SA10_gen.java
    #   (α,β-unsaturated carbonyls). Written here for cyclic esters rather than
    #   transcribed; the SA10 string itself is unused (see the disabled entry
    #   above).
    (
        "α,β-Unsaturated lactone (락톤 Michael acceptor)",
        "[#8;R]-,:[#6;R](=[OX1])-,:[#6;R;!$([#6]-[OX2H]);!$([#6]-[NX3])]"
        "=,:[#6;R;!$([#6]-[OX2H]);!$([#6]-[NX3])]",
    ),
    ("Aromatic nitroso (방향족 니트로소)", "[c][NX2]=[OX1]"),
    ("Alkyl azoxy (알킬 아조옥시)", "[#6][NX2]=[NX3+]([OX1-])[#6]"),
    # 2026-09-14 추가 — Benigni/Bossa 규칙베이스에서 보강.
    # 학습셋 결손 분석에서 지목되었고, 반응성 기전이 문헌상 명확하며
    # SMARTS 가 단순해 위양성 위험이 낮은 것만 골랐다.
    # 출처: Toxtree 3.1.0 소스 mutant/rules/ — 문자열 그대로 인용.
    #
    # SA22_gen — 결손 34건. 아자이드·트리아젠은 직접 DNA 알킬화제다.
    # 2026-09-29 수정: 두 이름이 SMARTS 와 서로 바뀌어 있었다.
    #   `[N]=[N]=[N]` 은 누적 이중결합이므로 아자이드 R-N=[N+]=[N-] 이고,
    #   `[N]=[N]-[N]` 은 트리아젠 R-N=N-NR2 다. 다카바진이 «Azide» 로,
    #   메틸 아자이드가 «Triazene» 으로 보고되고 있었다. SMARTS·발화 건수·
    #   정밀도는 그대로이며 표시 이름만 맞바꾼다 (성능 수치 영향 없음).
    # Attribution: verbatim from Toxtree 3.1.0 mutant/rules/SA22_gen.java.
    #   GPL-2.0, (c) Ideaconsult Ltd. Display names were swapped on 2026-09-29
    #   (see the note above); the SMARTS strings are unchanged from the source.
    ("Azide (아자이드)", "[N]=[N]=[N]"),
    ("Triazene (트리아젠)", "[N]=[N]-[N]"),
    # SA1_gen — 결손 12건. 아실 할라이드는 강한 친전자체다.
    # `[!$([OH1,SH1])]` 로 카복실산·티오산을 배제한다.
    # Attribution: verbatim from Toxtree 3.1.0 mutant/rules/SA1_gen.java.
    #   GPL-2.0, (c) Ideaconsult Ltd.
    ("Acyl halide (아실 할라이드)", "[!$([OH1,SH1])]C(=O)[Br,Cl,F,I]"),
)

# 리포트 Figure 2용 추가 패턴.
# 2026-09-03 수정: 구조 경보만으로는 ICH M7 Class 1/2를 확정할 수 없습니다
# (Class 1은 실측 Ames 양성 + 발암성 양성, Class 2는 실측 Ames 양성이 전제이며
# 경보만 있으면 Class 3입니다). 보고서 캡션에 그대로 나가던 등급 확정 표현을
# 제거하고 «structural alert»로 낮췄습니다. 방향족 아민·니트로·Michael acceptor는
# 위 Ashby-Tennant 목록과 SMARTS 문자열을 공유해 중복 계상되지 않습니다.
ICH_M7_REPORT_PATTERNS: tuple[tuple[str, str], ...] = (
    ("Aromatic amine (structural alert)", AROMATIC_AMINE_SMARTS),
    ("Nitro group (structural alert)", NITRO_AROMATIC_SMARTS),
    ("Michael acceptor (structural alert)", MICHAEL_ACCEPTOR_SMARTS),
    ("Heteroaromatic N (structural alert)", "[nH]1cccc1"),
    ("Anilide (aromatic amide)", "[c][NH]-[C]=[O]"),
    ("Ester (labile linkage)", "[#6][C](=O)O[#6]"),
)


# ---------------------------------------------------------------------------
# 알고리즘 경보 — SMARTS 로 표현할 수 없는 것
# ---------------------------------------------------------------------------
# 2026-09-14 추가: 다환방향족(PAH).
#
# Benigni/Bossa SA18_gen · SA19_gen 은 부분구조 매칭이 아니라 고리 융합 구조를
# 알고리즘으로 판정한다 (Toxtree 3.1.0 소스 mutant/rules/SA18_gen.java).
# 아래는 그 판정 논리를 RDKit 으로 옮긴 것이다.
#
#   1. 전체 고리 수가 3 미만이면 실패
#   2. 원자를 공유하는 고리끼리 묶는다 (융합 고리 집합)
#   3. 한 집합의 고리 수가 3 이상이어야 한다
#   4. 고리별 방향족성 판정
#        SA18 — (전체 원자 − 방향족 원자) < 2   관대
#        SA19 — (전체 원자 − 방향족 원자) == 0  엄격
#   5. 채택
#        SA18 — 헤테로방향족 고리 0개 AND 방향족 고리 > 2
#        SA19 — 헤테로방향족 고리 > 0  AND 방향족 고리 > 2
#
# 원 구현에 «더 견고한 방식(일반화 Hückel 규칙)으로 재구현할 것»이라는 TODO 가
# 달려 있다. 같은 방식을 쓰는 이상 그 한계도 함께 물려받는다.

PAH_MIN_FUSED_RINGS = 3


# ---------------------------------------------------------------------------
# 구조 경보 예외 (2026-09-29 신설)
# ---------------------------------------------------------------------------
# 경보는 "이 부분구조가 있으면 위험"이라는 단순 매칭이라 예외를 모릅니다.
# Derek Nexus 가 133종으로 하는 일의 상당 부분이 하나의 독성단을 치환 패턴별로
# 나눠 예외를 처리하는 것이고, 특이도 차이(79.9% 대 68.4%)가 거기서 옵니다.
#
# 아래 예외는 **학습셋(5,309건) 에서 도출한 뒤 기전을 확인해 고정한 것**입니다.
# 접히지 않은 Morgan 환경(반경 2)으로 후보를 뽑고, 지지 15건 이상 · 효과 20%p
# 이상 · Fisher 정확검정 BH 보정 q<0.05 · 골격 5-fold 중 4폴드 이상 동일 방향을
# 통과한 것만 남긴 뒤, 한계 FP:TP ≥ 4:1 로 선택하고 SMARTS 로 옮겼습니다.
# 검증셋(1,288건)은 최종 수치 보고에만 한 번 사용했습니다.
#
# 근거: 논문자료/02_검증결과/방향족아민_예외규칙_검토요청_20260929.md
#
# 각 항목: (예외 이름, 적용 SMARTS, 무효화 SMARTS 목록, 기전 메모)
#   적용 SMARTS 가 매칭되면 그 경보는 발화하지 않는다.
#   단, 무효화 SMARTS 중 하나라도 매칭되면 예외를 적용하지 않는다.
ALERT_EXCEPTIONS: dict[str, tuple[tuple[str, str, tuple[str, ...], str], ...]] = {
    "Aromatic amine (방향족 아민)": (
        (
            "para-sulfonyl aniline",
            "[NX3;H2,H1;!$(N[C,S,P]=[O,S]);!$(N-[N,O])]-c1ccc([SX4](=O)(=O))cc1",
            (
                # 아조 환원으로 다른 방향족 아민이 유리되므로 예외를 적용하지 않는다.
                "[c][NX2]=[NX2][c]",
            ),
            "강한 전자끌개인 설포닐기가 질소의 전자 밀도를 낮춰 P450 매개 "
            "N-수산화(나이트레늄 이온 경로의 필수 첫 단계)를 불리하게 만든다. "
            "설폰산은 생리적 pH 에서 음이온으로 존재해 균체 내 유입도 제한된다. "
            "학습셋 19건 전부 Ames 음성.",
        ),
        (
            "anthranilate ester",
            "[NX3;H2,H1;!$(N[C,S,P]=[O,S]);!$(N-[N,O])]-c1ccccc1[CX3](=O)[OX2][#6]",
            (
                "[c][NX2]=[NX2][c]",
                # 벤지딘·디아미노디페닐메탄형은 연결된 이중 아닐린 골격 자체가
                # 변이원성의 원인이므로 예외 대상에서 뺀다.
                "[NX3;H2,H1;!$(N[C,S,P]=[O,S]);!$(N-[N,O])]-c.[NX3;H2,H1;"
                "!$(N[C,S,P]=[O,S]);!$(N-[N,O])]-c",
            ),
            "오쏘 치환이 나이트레늄 이온 안정화에 필요한 평면 배치를 입체적으로 "
            "방해하고, 아미노기와 오쏘 카보닐 사이의 분자내 수소결합이 질소의 "
            "친핵성과 N-산화 접근성을 함께 떨어뜨린다. 카복실산·케톤까지 넓히면 "
            "안트라퀴논 아민이 딸려 들어오므로 에스터로 한정. 학습셋 13건 전부 "
            "Ames 음성.",
        ),
    ),
}


# 같은 독성단을 다른 이름으로 부르는 곳이 있어 별칭으로 묶는다.
# p66_ich_m7.CORE_M7_ALERT_SMARTS 는 폴백 목록이며 "Aromatic amine" 이라는
# 짧은 이름으로 **동일한 SMARTS** 를 한 번 더 들고 있다. 별칭이 없으면
# 예외로 경보를 끈 화합물이 그 폴백에서 다시 잡혀 예외가 무력화된다.
# (2026-09-29 — 실제로 그렇게 되어 s5 수치가 그대로였다.)
ALERT_NAME_ALIASES: dict[str, str] = {
    "Aromatic amine": "Aromatic amine (방향족 아민)",
}


def alert_exception_applies(mol: Chem.Mol, alert_name: str) -> str:
    """공개 진입점 — 다른 모듈의 경보 목록에서도 같은 예외를 적용한다."""
    return _exception_applies(mol, alert_name)


def _exception_applies(mol: Chem.Mol, alert_name: str) -> str:
    """이 경보에 대해 예외가 성립하면 예외 이름을, 아니면 빈 문자열을 반환."""
    rules = ALERT_EXCEPTIONS.get(alert_name)
    if not rules:
        rules = ALERT_EXCEPTIONS.get(ALERT_NAME_ALIASES.get(alert_name, ""))
    if not rules:
        return ""
    for exc_name, smarts, vetoes, _memo in rules:
        try:
            pat = Chem.MolFromSmarts(smarts)
            if pat is None or not mol.HasSubstructMatch(pat):
                continue
            blocked = False
            for v in vetoes:
                vpat = Chem.MolFromSmarts(v)
                if vpat is not None and mol.HasSubstructMatch(vpat):
                    blocked = True
                    break
            if not blocked:
                return exc_name
        except Exception:  # noqa: BLE001
            continue
    return ""


def _fused_ring_systems(mol: Chem.Mol) -> list[list[tuple[int, ...]]]:
    """원자를 공유하는 고리끼리 묶어 융합 고리 집합 목록을 돌려준다."""
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) < PAH_MIN_FUSED_RINGS:
        return []

    systems: list[list[tuple[int, ...]]] = []
    for ring in rings:
        atoms = set(ring)
        merged: list[tuple[int, ...]] = [ring]
        rest: list[list[tuple[int, ...]]] = []
        for system in systems:
            if any(atoms & set(r) for r in system):
                merged.extend(system)
            else:
                rest.append(system)
        systems = rest + [merged]
    return systems


def _ring_is_aromatic(mol: Chem.Mol, ring: tuple[int, ...], *, strict: bool) -> bool:
    """SA18(관대) / SA19(엄격) 의 고리 방향족성 판정."""
    non_aromatic = sum(
        1 for idx in ring if not mol.GetAtomWithIdx(idx).GetIsAromatic()
    )
    return non_aromatic == 0 if strict else non_aromatic < 2


def _ring_is_heterocyclic(mol: Chem.Mol, ring: tuple[int, ...]) -> bool:
    return any(mol.GetAtomWithIdx(idx).GetSymbol() != "C" for idx in ring)


def _scan_pah(mol: Chem.Mol) -> list[tuple[str, set[int]]]:
    """PAH 경보. (경보 이름, 하이라이트 원자) 목록을 돌려준다.

    Attribution: independent reimplementation of the decision logic of
    Toxtree 3.1.0 mutant/rules SA18 (carbocyclic polycyclic aromatics) and
    SA19 (heterocyclic polycyclic aromatics), written with RDKit ring
    perception. No Toxtree code or pattern was copied — these two alerts are
    the "2 algorithmic" patterns referred to in the article, and they are
    algorithmic precisely because the fused-ring and aromaticity conditions
    cannot be expressed as a single SMARTS.
    """
    hits: list[tuple[str, set[int]]] = []
    carbocyclic_atoms: set[int] = set()
    heterocyclic_atoms: set[int] = set()

    for system in _fused_ring_systems(mol):
        if len(system) < PAH_MIN_FUSED_RINGS:
            continue

        # SA18 — 관대한 방향족성, 헤테로 고리가 하나도 없어야 한다
        aromatic = [r for r in system if _ring_is_aromatic(mol, r, strict=False)]
        hetero = [r for r in aromatic if _ring_is_heterocyclic(mol, r)]
        if len(aromatic) > 2 and not hetero:
            for ring in aromatic:
                carbocyclic_atoms.update(ring)
            continue

        # SA19 — 엄격한 방향족성, 헤테로 고리가 하나 이상
        aromatic_strict = [r for r in system if _ring_is_aromatic(mol, r, strict=True)]
        hetero_strict = [r for r in aromatic_strict if _ring_is_heterocyclic(mol, r)]
        if len(aromatic_strict) > 2 and hetero_strict:
            for ring in aromatic_strict:
                heterocyclic_atoms.update(ring)

    if carbocyclic_atoms:
        hits.append((PAH_CARBOCYCLIC_NAME, carbocyclic_atoms))
    if heterocyclic_atoms:
        hits.append((PAH_HETEROCYCLIC_NAME, heterocyclic_atoms))
    return hits


PAH_CARBOCYCLIC_NAME = "Polycyclic aromatic hydrocarbon (다환방향족, 3고리 이상)"
PAH_HETEROCYCLIC_NAME = "Heterocyclic PAH (헤테로 다환방향족, 3고리 이상)"


def _scan_patterns(
    mol: Chem.Mol,
    patterns: tuple[tuple[str, str], ...],
    *,
    alert_names: list[str],
    highlight_atoms: set[int],
    seen_smarts: set[str],
) -> None:
    """패턴 목록을 훑어 매칭 이름·원자를 누적.

    2026-09-03 수정: 중복 제거 기준을 «이름»에서 «SMARTS 문자열»로 변경.
    이름만 다른 동일 SMARTS(예: 니트로기 vs ICH M7 Class 2)가 하나의 작용기를
    두세 번 계상해 알러트 개수를 부풀리던 문제.
    """
    for name, smarts in patterns:
        if smarts in seen_smarts:
            continue
        seen_smarts.add(smarts)
        try:
            pat = Chem.MolFromSmarts(smarts)
            if pat is None:
                continue
            matches = mol.GetSubstructMatches(pat)
            if matches:
                # 2026-09-29: 큐레이션된 예외에 해당하면 이 경보는 발화하지 않는다.
                if _exception_applies(mol, name):
                    continue
                if name not in alert_names:
                    alert_names.append(name)
                for match in matches:
                    highlight_atoms.update(match)
        except Exception:
            continue


def scan_ashby_tennant_alerts(mol: Chem.Mol | None) -> tuple[list[str], list[int]]:
    """
    매칭된 알러트 이름 목록과 하이라이트할 원자 인덱스 반환.
    process_reaction_engine alert_fn 호환.
    """
    if mol is None:
        return [], []

    alert_names: list[str] = []
    highlight_atoms: set[int] = set()

    _scan_patterns(
        mol,
        ASHBY_TENNANT_ALERTS,
        alert_names=alert_names,
        highlight_atoms=highlight_atoms,
        seen_smarts=set(),
    )

    # 2026-09-14 추가: SMARTS 로 표현할 수 없는 알고리즘 경보 (PAH).
    for name, atoms in _scan_pah(mol):
        if name not in alert_names:
            alert_names.append(name)
        highlight_atoms.update(atoms)

    return alert_names, sorted(highlight_atoms)


def _bonds_between_atoms(mol: Chem.Mol, atoms: list[int]) -> list[int]:
    bonds: set[int] = set()
    for i, atom_a in enumerate(atoms):
        for atom_b in atoms[i + 1 :]:
            bond = mol.GetBondBetweenAtoms(atom_a, atom_b)
            if bond is not None:
                bonds.add(bond.GetIdx())
    return sorted(bonds)


def scan_report_structural_alerts(mol: Chem.Mol | None) -> tuple[list[str], list[int], list[int]]:
    """
    리포트 Figure 2용 — Ashby-Tennant + ICH M7 SMARTS 통합 스캔.
    반환: (alert_names, highlight_atoms, highlight_bonds)
    """
    if mol is None:
        return [], [], []

    alert_names: list[str] = []
    atoms_set: set[int] = set()
    seen_smarts: set[str] = set()

    _scan_patterns(
        mol,
        ASHBY_TENNANT_ALERTS,
        alert_names=alert_names,
        highlight_atoms=atoms_set,
        seen_smarts=seen_smarts,
    )
    _scan_patterns(
        mol,
        ICH_M7_REPORT_PATTERNS,
        alert_names=alert_names,
        highlight_atoms=atoms_set,
        seen_smarts=seen_smarts,
    )

    atoms = sorted(atoms_set)
    return alert_names, atoms, _bonds_between_atoms(mol, atoms)


def get_alert_highlight_indices(
    mol: Chem.Mol | None,
    alert_names: list[str] | None = None,
) -> tuple[list[int], list[int]]:
    """
    rdMolDraw2D 하이라이트용 원자·결합 인덱스 반환.
    alert_names가 없거나 매칭 실패 시 통합 리포트 스캔으로 fallback.
    """
    if mol is None:
        return [], []

    names = list(alert_names or [])
    if names:
        # 2026-09-03 수정: 리포트 전용 패턴 이름도 조회 대상에 포함 (기존에는
        # Ashby-Tennant 이름만 있어 리포트 알러트가 항상 fallback으로 빠졌습니다).
        name_to_smarts = {**dict(ASHBY_TENNANT_ALERTS), **dict(ICH_M7_REPORT_PATTERNS)}
        highlight_atoms: set[int] = set()
        highlight_bonds: set[int] = set()

        for name in names:
            smarts = name_to_smarts.get(name)
            if not smarts:
                continue
            try:
                alert_query = Chem.MolFromSmarts(smarts)
                if alert_query is None:
                    continue
                for match in mol.GetSubstructMatches(alert_query):
                    atoms = list(match)
                    highlight_atoms.update(atoms)
                    for i, atom_a in enumerate(atoms):
                        for atom_b in atoms[i + 1 :]:
                            bond = mol.GetBondBetweenAtoms(atom_a, atom_b)
                            if bond is not None:
                                highlight_bonds.add(bond.GetIdx())
            except Exception:
                continue

        if highlight_atoms:
            return sorted(highlight_atoms), sorted(highlight_bonds)

    _, atoms, bonds = scan_report_structural_alerts(mol)
    return atoms, bonds


def alert_status_label(alerts: list[str]) -> tuple[str, str]:
    """(status_emoji, status_text) — POSITIVE / NEGATIVE."""
    if alerts:
        return "🔴", "POSITIVE — Structural alert(s) detected"
    return "🟢", "NEGATIVE — No Ashby-Tennant alert detected"
