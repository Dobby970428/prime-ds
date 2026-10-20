"""부형제-원료의약품 반응 규칙 — 확장판.

기존 `p72.AdvancedExcipientMatrixEngine`에는 SMIRKS 규칙이 두 개뿐입니다.

    Maillard_Imine_Bond   [CH1:1](=O).[NX3;H2,H1:2]>>[CH1:1]=[NX3:2]
    Formylation_Artifact  C(=O)O.[NX3;H2,H1:1]>>[NX3:1]C=O

부형제를 아무리 늘려도 규칙이 없으면 아무것도 안 잡힙니다. 이 모듈은
전처방 배합적합성에서 실제로 문제가 되는 반응을 추가합니다.

수록 근거는 Handbook of Pharmaceutical Excipients의 Incompatibilities
항목과 배합적합성 문헌에서 반복적으로 보고되는 경로를 기준으로 했습니다.

⚠️ SMARTS 대소문자 주의
    대문자 C = 지방족 탄소, 소문자 c = 방향족 탄소, [#6] = 모든 탄소
    지난번 분해 엔진에서 이 때문에 아릴 에스터가 전부 누락됐습니다.
    여기서는 방향족을 포함해야 하는 자리에 [#6]을 씁니다.

2026-08-25 작성.
"""

from __future__ import annotations

from typing import Any

__all__ = [
    "EXTENDED_RULES",
    "RULE_INFO",
    "build_reaction_rules",
    "rule_label",
]


# ---------------------------------------------------------------------------
# 규칙 정의
#
# SMIRKS는 (부형제 조각, 원료의약품) 순서로 적용됩니다. 엔진이 순서를
# 바꿔 한 번 더 시도하므로 방향은 크게 중요하지 않지만, 읽기 쉽게
# 부형제를 앞에 두었습니다.
# ---------------------------------------------------------------------------

EXTENDED_RULES: tuple[tuple[str, str], ...] = (
    # ── Maillard ─────────────────────────────────────────────────────────
    #
    # ⚠️ 기존 규칙은 `[NX3;H2,H1]`만 요구해서 두 가지가 잘못됐습니다.
    #
    #   ① 1급 **아마이드** 질소에도 이민을 붙였습니다.
    #      아테놀올에서 `CC(=O)N=CC(O)CO`가 나오던 것이 이 경우입니다.
    #      아마이드 질소는 친핵성이 약해 제형 조건에서 Schiff base를
    #      만들지 않습니다.
    #
    #   ② 2급 아민(H1)에도 중성 이민 `C=N`을 만들려 했습니다.
    #      2급 아민이 이민을 만들면 질소 원자가가 초과됩니다.
    #      2급 아민은 이미늄을 거쳐 **엔아민**이 됩니다.
    #
    # 그래서 1급 아민 전용 이민 규칙과 2급 아민 전용 엔아민 규칙으로
    # 나눴습니다.
    #
    # 2026-09-03 수정: 산물 쪽 H 개수 고정을 전 규칙에서 제거했습니다.
    #   Michael·과산화물 규칙에서 실제로 라디칼을 만든 것과 **같은 종류의
    #   위험**이라, 아직 드러나지 않은 규칙도 함께 풀었습니다. 대표 기질로
    #   전후 산물이 동일함을 확인했습니다(H 재계산 결과가 같음).
    (
        "Maillard_Imine_Bond",
        "[CX3H1:1](=[OX1])."
        "[NX3;H2;!$(NC=O);!$(N=*);!$(N[N,O,S]);!$(Nc);!$(NS(=O)=O):2]"
        ">>[C:1]=[NX2:2]",
    ),
    (
        "Maillard_Enamine",
        "[CX3H1:1](=[OX1])[CX4;H1,H2,H3:2]."
        "[NX3;H1;!$(NC=O);!$(N=*);!$(N[N,O,S]);!$(Nc);!$(NS(=O)=O):3]"
        ">>[NX3:3][C:1]=[C:2]",
    ),
    (
        "Formylation_Artifact",
        "[CX3H1:1](=[OX1])[OX2H1]."
        "[NX3;H2,H1;!$(NC=O);!$(N=*);!$(N[N,O,S]):2]>>[NX3:2][C:1]=O",
    ),

    # ── 아마이드 생성 ────────────────────────────────────────────────────
    # 카복실산 부형제(스테아르산·구연산·타르타르산 등) + 1·2급 아민 API
    # 고체상에서도 고온·장기보관 시 보고되는 경로입니다.
    (
        "Amide_Condensation",
        "[CX3:1](=[OX1:2])[OX2H1].[NX3;H2,H1;!$(NC=O);!$(N=*);!$(N[N,O]):3]"
        ">>[CX3:1](=[OX1:2])[NX3:3]",
    ),

    # ── 에스터 생성 · 에스터 교환 ────────────────────────────────────────
    # 카복실산 부형제 + 알코올/페놀 API
    (
        "Esterification",
        "[CX3:1](=[OX1:2])[OX2H1].[OX2H1;$([OX2H1][#6]):3]>>[CX3:1](=[OX1:2])[OX2:3]",
    ),
    # 에스터 부형제(스테아레이트·파라벤·폴리소르베이트) + 알코올 API
    (
        "Transesterification",
        "[CX3:1](=[OX1:2])[OX2][#6].[OX2H1;$([OX2H1][#6]):3]>>[CX3:1](=[OX1:2])[OX2:3]",
    ),
    # 에스터 부형제 + 1·2급 아민 API → 아마이드 (aminolysis)
    (
        "Ester_Aminolysis",
        "[CX3:1](=[OX1:2])[OX2][#6].[NX3;H2,H1;!$(NC=O);!$(N=*):3]"
        ">>[CX3:1](=[OX1:2])[NX3:3]",
    ),

    # ── Michael 부가 ─────────────────────────────────────────────────────
    # α,β-불포화 부형제(푸마르산·소르빈산·말레산) + 아민 API
    #
    # 2026-09-03 수정: 산물의 H 개수 고정(`[CX3H1:2]`)을 없애고 `[C:2]`로 바꿨습니다.
    #   부가가 끝나면 α·β 탄소는 둘 다 sp3가 되는데, 템플릿이 β-탄소의 H를
    #   1개로 못박아 원자가가 하나 모자란 **탄소 라디칼**(rad=1)이 나왔습니다.
    #   H 지정을 빼면 RDKit이 결합 상태를 보고 다시 계산합니다.
    (
        "Michael_Addition_Amine",
        "[CX3:1]=[CX3:2][CX3:3]=[OX1:4].[NX3;H2,H1;!$(NC=O);!$(N=*):5]"
        ">>[C:1]([NX3:5])[C:2][CX3:3]=[OX1:4]",
    ),
    # 티올 API (캅토프릴 등)
    (
        "Michael_Addition_Thiol",
        "[CX3:1]=[CX3:2][CX3:3]=[OX1:4].[SX2H1:5]"
        ">>[C:1]([SX2:5])[C:2][CX3:3]=[OX1:4]",
    ),

    # ── 중아황산 부가 ────────────────────────────────────────────────────
    # 메타중아황산나트륨·아황산염 + 알데하이드/케톤 API
    # 항산화제로 넣은 아황산염이 카보닐 API와 부가물을 만드는 잘 알려진 경로.
    # 2026-09-03 수정: 산물 H 고정 제거 (`[CX4H1:5]`→`[C:5]`, `[OX2H1:6]`→`[O:6]`).
    (
        "Bisulfite_Adduct",
        "[OX2H1:3][SX3:1](=[OX1:2])[OX2H1].[CX3H1:5]=[OX1:6]"
        ">>[C:5]([O:6])[SX4:1](=[OX1:2])(=[OX1])[O:3]",
    ),
    (
        "Bisulfite_Adduct_Anion",
        "[OX1-:3][SX3:1](=[OX1:2])[OX2H1].[CX3H1:5]=[OX1:6]"
        ">>[C:5]([O:6])[SX4:1](=[OX1:2])(=[OX1])[O-:3]",
    ),

    # ── 과산화물 산화 ────────────────────────────────────────────────────
    # 포비돈·크로스포비돈·PEG·폴리소르베이트의 잔류 과산화물
    # 3급 아민 → N-oxide
    #
    # 2026-09-03 수정: 두 번째 산물을 `[OX2H1:1]` → `[O:1]`(물)로 바꿨습니다.
    #   H2O2 + R3N → R3N⁺-O⁻ + H2O 이므로 부산물은 물인데, `[OX2H1:1]`이
    #   H를 1개로 못박아 **수산화 라디칼 [OH]**(rad=1)가 나왔습니다.
    #   H 지정을 빼면 RDKit이 H를 2개로 채워 물이 됩니다.
    (
        "Peroxide_N_Oxidation",
        "[OX2H1:1][OX2H1].[NX3;H0;$([NX3]([#6])([#6])[#6]);!$(N[C,S]=[O,S,N]);!$(N=*):2]"
        ">>[N+:2][O-].[O:1]",
    ),
    # 티오에터 → 설폭사이드
    (
        "Peroxide_S_Oxidation",
        "[OX2H1:1][OX2H1].[SX2;$([SX2]([#6])[#6]);!$(S=*):2]>>[S+:2][O-].[O:1]",
    ),

    # ── 아세탈 · 헤미아세탈 ──────────────────────────────────────────────
    # 알데하이드 부형제 잔류물 + 알코올 API
    # 2026-09-03 수정: 산물 H 고정 제거 (`[CX4H1:1]`→`[C:1]`, `[OX2H1:2]`→`[O:2]`).
    (
        "Hemiacetal_Formation",
        "[CX3H1:1]=[OX1:2].[OX2H1;$([OX2H1][#6]):3]>>[C:1]([O:2])[OX2:3]",
    ),
)


# ---------------------------------------------------------------------------
# 규칙 설명 — 화면·보고서 표시용
# ---------------------------------------------------------------------------

RULE_INFO: dict[str, dict[str, str]] = {
    "Maillard_Imine_Bond": {
        "en": "Maillard imine (aldehyde + amine)",
        "ko": "Maillard 이민 (알데하이드 + 아민)",
        "group": "maillard",
        "note_en": "Reducing sugars condense with a primary amine to give a Schiff base. "
                   "Amide nitrogens are excluded — they are poor nucleophiles and do not form "
                   "imines under formulation conditions.",
        "note_ko": "환원당이 1급 아민과 축합해 Schiff base를 만듭니다. 아마이드 질소는 "
                   "친핵성이 약해 제형 조건에서 이민을 만들지 않으므로 제외했습니다.",
    },
    "Maillard_Enamine": {
        "en": "Maillard enamine (aldehyde + secondary amine)",
        "ko": "Maillard 엔아민 (알데하이드 + 2급 아민)",
        "group": "maillard",
        "note_en": "A secondary amine cannot form a neutral imine — it proceeds through the "
                   "iminium to the enamine. Requires an α-hydrogen on the aldehyde.",
        "note_ko": "2급 아민은 중성 이민을 만들 수 없고 이미늄을 거쳐 엔아민이 됩니다. "
                   "알데하이드에 α-수소가 있어야 합니다.",
    },
    "Formylation_Artifact": {
        "en": "Formylation (formic acid + amine)",
        "ko": "포르밀화 (포름산 + 아민)",
        "group": "formylation",
        "note_en": "Formic acid, typically from macrogol or polysorbate oxidation, formylates "
                   "primary/secondary amines.",
        "note_ko": "마크로골·폴리소르베이트 산화로 생긴 포름산이 1·2급 아민을 포르밀화합니다.",
    },
    "Amide_Condensation": {
        "en": "Amide condensation (carboxylic acid + amine)",
        "ko": "아마이드 축합 (카복실산 + 아민)",
        "group": "acylation",
        "note_en": "Acidic excipients such as stearic, citric or tartaric acid can acylate an "
                   "amine drug substance on prolonged storage.",
        "note_ko": "스테아르산·구연산·타르타르산 같은 산성 부형제가 장기 보관 중 아민 원료를 아실화합니다.",
    },
    "Esterification": {
        "en": "Esterification (carboxylic acid + alcohol)",
        "ko": "에스터화 (카복실산 + 알코올)",
        "group": "acylation",
        "note_en": "Acidic excipient esterifies a hydroxyl group on the drug substance.",
        "note_ko": "산성 부형제가 원료의 수산기를 에스터화합니다.",
    },
    "Transesterification": {
        "en": "Transesterification (ester excipient + alcohol)",
        "ko": "에스터 교환 (에스터 부형제 + 알코올)",
        "group": "acylation",
        "note_en": "Ester excipients such as stearates, parabens or polysorbates exchange with "
                   "a hydroxyl group on the drug substance.",
        "note_ko": "스테아레이트·파라벤·폴리소르베이트 같은 에스터 부형제가 원료의 수산기와 교환됩니다.",
    },
    "Ester_Aminolysis": {
        "en": "Ester aminolysis (ester excipient + amine)",
        "ko": "에스터 아미놀리시스 (에스터 부형제 + 아민)",
        "group": "acylation",
        "note_en": "An amine drug substance attacks an ester excipient to give the corresponding "
                   "amide, liberating the alcohol.",
        "note_ko": "아민 원료가 에스터 부형제를 공격해 아마이드를 만들고 알코올을 유리합니다.",
    },
    "Michael_Addition_Amine": {
        "en": "Michael addition (α,β-unsaturated + amine)",
        "ko": "Michael 부가 (α,β-불포화 + 아민)",
        "group": "michael",
        "note_en": "Fumarate, maleate and sorbate excipients act as Michael acceptors toward "
                   "amine drug substances.",
        "note_ko": "푸마르산·말레산·소르빈산 부형제가 아민 원료에 대해 Michael 수용체로 작용합니다.",
    },
    "Michael_Addition_Thiol": {
        "en": "Michael addition (α,β-unsaturated + thiol)",
        "ko": "Michael 부가 (α,β-불포화 + 티올)",
        "group": "michael",
        "note_en": "Thiol drug substances add rapidly to unsaturated excipients.",
        "note_ko": "티올 원료는 불포화 부형제에 빠르게 부가됩니다.",
    },
    "Bisulfite_Adduct": {
        "en": "Bisulfite adduct (sulfite + carbonyl)",
        "ko": "중아황산 부가 (아황산염 + 카보닐)",
        "group": "bisulfite",
        "note_en": "Sulfite antioxidants form α-hydroxysulfonate adducts with aldehyde and "
                   "ketone drug substances — a well documented incompatibility.",
        "note_ko": "항산화제로 쓰는 아황산염이 알데하이드·케톤 원료와 α-하이드록시설포네이트 "
                   "부가물을 만듭니다. 문헌에 잘 정리된 배합금기입니다.",
    },
    "Bisulfite_Adduct_Anion": {
        "en": "Bisulfite adduct (bisulfite anion + carbonyl)",
        "ko": "중아황산 부가 (중아황산 이온 + 카보닐)",
        "group": "bisulfite",
        "note_en": "Same pathway as the neutral form, written for the bisulfite anion.",
        "note_ko": "중성형과 같은 경로이며 중아황산 이온 표기를 처리합니다.",
    },
    "Peroxide_N_Oxidation": {
        "en": "Peroxide N-oxidation (tertiary amine → N-oxide)",
        "ko": "과산화물 N-산화 (3급 아민 → N-oxide)",
        "group": "oxidation",
        "note_en": "Residual peroxide in povidone, macrogol or polysorbate oxidises a tertiary "
                   "amine to the N-oxide.",
        "note_ko": "포비돈·마크로골·폴리소르베이트의 잔류 과산화물이 3급 아민을 N-oxide로 산화합니다.",
    },
    "Peroxide_S_Oxidation": {
        "en": "Peroxide S-oxidation (thioether → sulfoxide)",
        "ko": "과산화물 S-산화 (티오에터 → 설폭사이드)",
        "group": "oxidation",
        "note_en": "Residual peroxide oxidises a thioether to the sulfoxide.",
        "note_ko": "잔류 과산화물이 티오에터를 설폭사이드로 산화합니다.",
    },
    "Hemiacetal_Formation": {
        "en": "Hemiacetal formation (aldehyde + alcohol)",
        "ko": "헤미아세탈 생성 (알데하이드 + 알코올)",
        "group": "maillard",
        "note_en": "Residual aldehyde adds to a hydroxyl group on the drug substance.",
        "note_ko": "잔류 알데하이드가 원료의 수산기에 부가됩니다.",
    },
}


def rule_label(name: str, lang: str = "en") -> str:
    info = RULE_INFO.get(name)
    if not info:
        return name.replace("_", " ")
    return info.get("ko" if lang == "ko" else "en", name)


def rule_group(name: str) -> str:
    return (RULE_INFO.get(name) or {}).get("group", "other")


def build_reaction_rules() -> dict[str, Any]:
    """SMIRKS를 RDKit 반응 객체로. 파싱 실패한 규칙은 조용히 빠집니다."""
    from rdkit.Chem import AllChem

    rules: dict[str, Any] = {}
    for name, smirks in EXTENDED_RULES:
        try:
            rxn = AllChem.ReactionFromSmarts(smirks)
        except Exception:
            continue
        if rxn is None:
            continue
        try:
            rxn.Initialize()
        except Exception:
            pass
        rules[name] = rxn
    return rules


# ---------------------------------------------------------------------------
# 검증용 대표 기질 — (부형제 조각, 원료의약품)
#
# 2026-09-03 추가. `validate()`가 파싱만 보고 **산물은 보지 않아서**
# Michael 부가·과산화물 산화가 라디칼 산물을 내놓는데도 bad=[] 였습니다.
# 규칙마다 대표 쌍을 두고 실제로 돌려 산물의 sanitize·라디칼을 검사합니다.
# ---------------------------------------------------------------------------

VALIDATION_SUBSTRATES: dict[str, tuple[str, str]] = {
    # 글리세르알데하이드 + 펜에틸아민 (1급 아민)
    "Maillard_Imine_Bond": ("O=CC(O)CO", "NCCc1ccccc1"),
    # 글리세르알데하이드 + N-메틸펜에틸아민 (2급 아민)
    "Maillard_Enamine": ("O=CC(O)CO", "CNCCc1ccccc1"),
    "Formylation_Artifact": ("OC=O", "NCCc1ccccc1"),
    "Amide_Condensation": ("CCCCCCCCCCCCCCCCCC(=O)O", "NCCc1ccccc1"),
    "Esterification": ("OC(=O)CC(O)(CC(=O)O)C(=O)O", "OCCc1ccccc1"),
    "Transesterification": ("COC(=O)c1ccc(O)cc1", "OCCc1ccccc1"),
    "Ester_Aminolysis": ("COC(=O)c1ccc(O)cc1", "NCCc1ccccc1"),
    # 푸마르산 + N-메틸펜에틸아민
    "Michael_Addition_Amine": ("OC(=O)/C=C/C(=O)O", "CNCCc1ccccc1"),
    # 푸마르산 + 캅토프릴 (티올)
    "Michael_Addition_Thiol": ("OC(=O)/C=C/C(=O)O", "CC(CS)C(=O)N1CCCC1C(=O)O"),
    "Bisulfite_Adduct": ("OS(=O)O", "O=Cc1ccccc1"),
    "Bisulfite_Adduct_Anion": ("[O-]S(=O)O", "O=Cc1ccccc1"),
    # 과산화수소 + 3급 아민(N,N-디메틸펜에틸아민)
    "Peroxide_N_Oxidation": ("OO", "CN(C)CCc1ccccc1"),
    # 과산화수소 + 티오에터
    "Peroxide_S_Oxidation": ("OO", "CSCCc1ccccc1"),
    "Hemiacetal_Formation": ("O=Cc1ccccc1", "OCCc1ccccc1"),
}


def _inspect_products(rxn: Any, pair: tuple[str, str]) -> tuple[list[str], list[str]]:
    """대표 쌍을 실제로 돌려 (산물 SMILES 목록, 문제 목록)을 냅니다."""
    from rdkit import Chem

    reactant_a = Chem.MolFromSmiles(pair[0])
    reactant_b = Chem.MolFromSmiles(pair[1])
    if reactant_a is None or reactant_b is None:
        return [], ["validation substrate failed to parse"]

    try:
        product_sets = rxn.RunReactants((reactant_a, reactant_b))
    except Exception as exc:
        return [], [f"RunReactants raised {type(exc).__name__}: {exc}"]
    if not product_sets:
        return [], ["no product from representative substrates"]

    smiles_out: list[str] = []
    problems: list[str] = []
    for prod_tuple in product_sets:
        for prod in prod_tuple:
            try:
                Chem.SanitizeMol(prod)
            except Exception as exc:
                problems.append(f"sanitize failed: {type(exc).__name__}: {exc}")
                continue
            radicals = sum(a.GetNumRadicalElectrons() for a in prod.GetAtoms())
            smi = Chem.MolToSmiles(prod)
            if radicals:
                problems.append(f"radical product ({radicals} e⁻): {smi}")
            if Chem.MolFromSmiles(smi) is None:
                problems.append(f"product SMILES does not round-trip: {smi}")
            if smi not in smiles_out:
                smiles_out.append(smi)
    return smiles_out, problems


def validate() -> dict[str, Any]:
    """규칙 파싱 + **산물 유효성**까지 확인 — 검증용.

    2026-09-03 수정: 전에는 `ReactionFromSmarts`가 None을 내는지만 봤습니다.
    파싱은 되지만 산물 템플릿의 H 고정 때문에 원자가가 모자란 라디칼이
    나오는 규칙(Michael 부가 2종·과산화물 산화 2종)을 전혀 잡지 못했습니다.
    이제 대표 기질로 실제 반응을 돌려 sanitize·라디칼·round-trip을 봅니다.

    반환 키(`total`·`ok`·`bad`)는 하위 호환을 위해 유지하고,
    `product_checks`·`radical_free`를 덧붙였습니다.
    """
    from rdkit.Chem import AllChem

    ok, bad = [], []
    product_checks: dict[str, dict[str, Any]] = {}

    for name, smirks in EXTENDED_RULES:
        try:
            rxn = AllChem.ReactionFromSmarts(smirks)
        except Exception as exc:
            bad.append((name, str(exc)))
            continue
        if rxn is None:
            bad.append((name, "returned None"))
            continue

        pair = VALIDATION_SUBSTRATES.get(name)
        if pair is None:
            ok.append(name)
            product_checks[name] = {"skipped": "no representative substrate"}
            continue

        try:
            rxn.Initialize()
        except Exception:
            pass

        products, problems = _inspect_products(rxn, pair)
        product_checks[name] = {
            "substrates": list(pair),
            "products": products,
            "problems": problems,
        }
        if problems:
            bad.append((name, "; ".join(problems)))
        else:
            ok.append(name)

    return {
        "total": len(EXTENDED_RULES),
        "ok": ok,
        "bad": bad,
        "product_checks": product_checks,
        "radical_free": not bad,
    }
