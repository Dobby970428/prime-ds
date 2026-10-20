"""
화학 반응식(SMIRKS) 기반 불순물 자동 예측 엔진.

RDKit rdChemReactions 모듈로 NDSRI(니트로사민) 형성 반응을 시뮬레이션하고,
DB SMILES 대조로 실측/이론 분류를 수행합니다.
"""

# 경로 A

# 니트로소화 담당(니트로사민+시약 분류)

# ① 시약이 뭔지 분류함(캔버스에 그린 시약 구조를 보고 역할 판정)
# 아질산염(nitrite)   →  니트로소화 가능
# 산화제              →  N-oxide 생성 가능
# 환원제              →  니트로기 환원 가능
# 알킬화제            →  알킬 부산물 가능
# 산                  →  가수분해 촉진

# ② 니트로사민을 만듦
# 2급 아민 + 아질산염  →  N-니트로사민
# 3급 아민 + 아질산염  →  탈알킬화 후 니트로사민

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from rdkit import Chem
from rdkit.Chem import Descriptors, rdChemReactions, rdMolDescriptors

# 아질산염(nitrite) 공정 조건 — NO2- / HONO 소스
NITRITE_SMILES = "[N+](=O)[O-]"

# ---------------------------------------------------------------------------
# [폐기 2026-08-26] 요청 SMIRKS — 실행 규칙에서 뺐습니다.
#
#   "[N:1]([C:2])[H:3].[N:4](=[O:5])[O-]>>[N:1]([C:2])[N:4]=[O:5]"
#
# 결함이 셋이었습니다.
#   ① 전하  — 아질산 이온의 [O-]만 떼고 질소의 +1을 그대로 물려받아
#             상대 음이온 없는 [N+]=O 가 생성됩니다. 성립하지 않는 구조입니다.
#   ② 수소  — [H:3] 매핑 때문에 AddHs 된 분자가 그대로 출력되어
#             DB의 표준 SMILES와 대조되지 않고 화면에서도 읽히지 않습니다.
#   ③ 아미드 — [N:1]([C:2])[H:3] 에 아미드 배제 조건이 없어
#             아테놀올의 1차 아마이드까지 니트로소화했습니다
#             (C(=O)N([H])[N+]=O).
#
# 같은 산물의 올바른 형태는 SECONDARY_AMINE_NITROSATION_SMIRKS 가
# 이미 정확히 만들어 냅니다. 즉 이 규칙은 중복이면서 형태만 잘못된
# 산물을 더하고 있었습니다. 제거해도 잃는 산물이 없음을 확인했습니다
# (원료 10종 × 부형제 8종 전수 대조, md파일/p30_니트로소화_규칙_수정.md).
#
# 기록을 위해 문자열만 남깁니다. 어디에서도 실행되지 않습니다.
# ---------------------------------------------------------------------------
RETIRED_SECONDARY_AMINE_NITROSATION_SMIRKS_USER = (
    "[N:1]([C:2])[H:3].[N:4](=[O:5])[O-]>>[N:1]([C:2])[N:4]=[O:5]"
)

# ---------------------------------------------------------------------------
# RDKit 경고 "reactant 1 has no mapped atoms." 에 대하여  (2026-09-15 기록)
#
#   아래 두 니트로소화 규칙은 두 번째 반응물 `[N+](=O)[O-]`(아질산염)에
#   원자 매핑을 두지 않는다. RDKit 은 반응물의 원자를 산물로 추적할 수
#   없을 때 이 경고를 내므로, 니트로소화가 시도될 때마다 한 번씩 찍힌다.
#   s9 실행에서 수십 번 나오는 것은 이 때문이다.
#
#   **버그가 아니라 의도된 형태다.** 산물 템플릿의 `N=O` 도 매핑이 없어,
#   RDKit 이 아질산염 원자를 옮기는 대신 산물 쪽에서 새로 만든다.
#   아질산염 원자는 버려진다. 결과 구조는 동일하다.
#
#   같은 기법을 TERTIARY_AMINE_NITROSATION_SMIRKS 가 이탈 알킬기에도
#   쓴다(아래 주석 참조) — 매핑을 빼서 원자를 지우는 SMIRKS 관용구다.
#
#   매핑을 넣은 판본은 RETIRED_..._USER 로 위에 남아 있으나, 아테놀올의
#   1차 아마이드까지 니트로소화하는 등의 문제로 은퇴시켰고 현재 형태를
#   원료 10종 × 부형제 8종 전수 대조로 검증했다
#   (md파일/p30_니트로소화_규칙_수정.md).
#
#   → 경고를 없애려고 매핑을 되살리지 말 것. 되살린다면 같은 전수 대조를
#     다시 수행해 산물이 동일함을 확인한 뒤에 할 것.
# ---------------------------------------------------------------------------

# [등록] 2급 아민(R2NH) 니트로사민화 — rdChemReactions 구동용
# 2026-09-03 수정: 반응물 질소를 [#7](방향족 n 포함) → [NX3](지방족 3배위 질소)로 좁혔습니다.
#   아래 SECONDARY_AMINE_SMARTS 주석대로 "방향족 고리 NH는 아민이 아니므로 기질이 아님"인데
#   SMIRKS만 [#7]이라 인돌·피롤 NH까지 니트로소화했습니다(수마트립탄에서 인돌 N-니트로소 위양성).
SECONDARY_AMINE_NITROSATION_SMIRKS = (
    "[NX3;H1;D2;!$([#7][#7])!$([#7]=*)!$([#7][C,S]=[!#6,#8,#16]):1]"
    "([#6:2])[#6:3].[N+](=O)[O-]>>[NX3;H0;D3:1]([#6:2])([#6:3])N=O"
)

# 2급 아민 작용기 탐지 (아미드·설폰아미드·히드라진·이민 제외)
#   NX3 = 지방족 질소만. 소문자 n(피롤·인돌 등 방향족 고리 질소)은 제외된다.
#   방향족 고리 NH는 아민이 아니므로 니트로사민화 기질이 아니다.
SECONDARY_AMINE_SMARTS = (
    "[NX3;H1;D2;!$(N[N,O,S]);!$(N=*);!$(N[C,S]=[O,N,S]);!$(NS(=O)=O)]"
    "([CX4,c])[CX4,c]"
)

# 3급 아민 (니트로소 탈알킬화로 NDSRI 생성 가능)
#   NX3 = 지방족 질소만. 카페인처럼 고리 질소가 D3/H0인 경우를 배제한다.
TERTIARY_AMINE_SMARTS = (
    "[NX3;H0;D3;!$(N[N,O,S]);!$(N=*);!$(N[C,S]=[O,N,S]);!$(NS(=O)=O);!$([N+])]"
    "([CX4,c])([CX4,c])[CX4,c]"
)

# 아질산염 / NOx 시약 탐지 (캔버스에 그려진 시약 구조)
NITRITE_SOURCE_SMARTS: tuple[str, ...] = (
    "[N+](=O)[O-]",
    "[N](=O)O",
    "O=N[OH]",
    "[N](=O)=O",
    "N#N",
    "[Na+].[N+](=O)[O-]",
)

# 시약 역할 분류 SMARTS
REAGENT_ROLE_PATTERNS: dict[str, tuple[str, ...]] = {
    "nitrite": NITRITE_SOURCE_SMARTS,
    "oxidant": (
        "[O][O]",
        "[O-][O+]",
        "[Mn+7]",
        "O=O",
        "[N+](=O)([O-])[O-]",
    ),
    "reductant": (
        "[BH3]",
        "[BH4-]",
        "[AlH]",
        "[SiH]",
    ),
    # 2026-09-03 수정: 단독 "[Br]"·"[I]"는 브로민화물 이온·아릴 브로마이드까지 전부 알킬화제로
    #   잡았습니다(클로로벤젠도 alkylating 판정). sp3 탄소에 결합한 이탈기 · 에폭사이드 ·
    #   설포네이트 에스터로 좁혔습니다.
    "alkylating": (
        "[CX4][Cl,Br,I]",
        "C1OC1",
        "[CX4]OS(=O)(=O)",
        "[CX4][OX2]S(=O)(=O)[c,C]",
    ),
    # 2026-09-03 수정: "[OX2H][CX3](=O)[OX2H]"는 카보닐 탄소에 OH가 둘 붙은 탄산에만 매치되어
    #   아세트산·벤조산이 False였습니다(사실상 죽은 조건). 일반 카복실산 패턴으로 교체.
    "acid": (
        "[CX3](=O)[OX2H1]",
        "S(=O)(=O)[OH]",
        "[OH][C,c]=O",
    ),
}

# 3급 아민 + 아질산염 → 니트로소 탈알킬화 → N-니트로소아민 (NDSRI)
#
#   R2N-CH2R'  +  HNO2  →  R2N-N=O  +  R'CHO
#
#   3급 아민은 질소에 수소가 없으므로 곧바로 니트로소화될 수 없다.
#   실제 경로는 니트로소암모늄 이온을 거쳐 알킬기 하나가 카르보닐로 떨어져 나가고
#   그 자리에 니트로소기가 들어가는 "니트로소 탈알킬화(nitrosative dealkylation)"다.
#   따라서 이탈하는 알킬기는 α-수소를 가져야 한다 → [CX4;H1,H2,H3]
#
#   이탈 알킬기는 생성물 템플릿에 매핑하지 않아 제거된다.
#   동반 생성되는 카르보닐(포름알데히드·아세트알데히드 등)은 별도 산물로 내보내지 않는다.
TERTIARY_AMINE_NITROSATION_SMIRKS = (
    "[NX3;H0;D3;!$(N[N,O,S]);!$(N=*);!$(N[C,S]=[O,N,S]);!$(NS(=O)=O);!$([N+]):1]"
    "([CX4,c:2])([CX4,c:3])[CX4;H1,H2,H3]"
    ".[N+](=O)[O-]"
    ">>[N:1]([#6:2])([#6:3])[N]=[O]"
)

# 시약 기반 부반응 SMIRKS (substrate + reagent)
REAGENT_SIDE_REACTIONS: tuple[tuple[str, str, str, str], ...] = (
    (
        "reagent_ndsri_secondary",
        "2급 아민 + 시약(아질산) → N-니트로소amine",
        "nitrite",
        SECONDARY_AMINE_NITROSATION_SMIRKS,
    ),
    (
        "reagent_ndsri_tertiary",
        "3급 아민 + 시약(아질산) → 니트로소 탈알킬화 → NDSRI",
        "nitrite",
        TERTIARY_AMINE_NITROSATION_SMIRKS,
    ),
    (
        "reagent_oxidation_aromatic",
        "방향족 N-산화 (산화제 시약)",
        "oxidant",
        # 2026-09-03 수정: [nH:1](피롤형 질소)은 N-oxide를 만들지 않고, 정작 N-oxide를
        #   만드는 피리딘형 질소(n;X2;H0)는 걸리지 않았습니다. 대상을 정반대로 잡던 것을 교정.
        "[n;X2;H0;!$([n+]):1]>>[n+:1][O-]",
    ),
    (
        "reagent_alkylation_amine",
        "2급 아민 알킬화 부반응 (알킬화 시약)",
        "alkylating",
        # 2026-09-03 수정: 산물의 [#6:5]가 반응물 어디에도 없어 RDKit이 새 메틸 탄소를 만들어
        #   붙였습니다(어떤 알킬화제를 넣어도 항상 메틸화). 반응물에 [CX4:5][Cl,Br,I]를 두고
        #   산물이 그 탄소를 참조하도록 교정 — 아릴 할라이드는 CX4가 아니므로 반응하지 않습니다.
        "[NX3;H1;D2;!$([#7][#7])!$([#7]=*);!$([#7][C,S]=[O,N,S]):1]([#6:2])[#6:3]"
        ".[CX4:5][Cl,Br,I]"
        ">>[NX3;H0;D3:1]([#6:2])([#6:3])[CX4:5]",
    ),
    (
        "reagent_reduction_carbonyl",
        "카르보닐 환원 부반응 (환원 시약)",
        "reductant",
        "[#6:1](=O)[#6:2]>>[#6:1](O)[#6:2]",
    ),
    (
        "reagent_ester_hydrolysis_acid",
        "에스터 가수분해 (산성 시약)",
        "acid",
        # 2026-09-03 수정: 산물을 카복실레이트 음이온([O-])으로 만들면 PubChem DB의 중성형
        #   canonical SMILES와 절대 일치하지 않아 항상 "이론적"으로 분류됩니다. 중성형으로 통일.
        "[#6:1](=O)O[#6:2]>>[#6:1](=O)[OX2H].[#6:2][OH]",
    ),
)

REACTION_REGISTRY: tuple[tuple[str, str, str], ...] = (
    (
        "secondary_amine_nitrosation",
        "2급 아민 + 아질산염 → N-니트로소amine (NDSRI)",
        SECONDARY_AMINE_NITROSATION_SMIRKS,
    ),
)


class ImpurityClassification(str, Enum):
    VALIDATED = "validated"
    THEORETICAL = "theoretical"


@dataclass(frozen=True)
class AmineSiteInfo:
    secondary_count: int
    tertiary_count: int = 0

    @property
    def has_secondary_amine(self) -> bool:
        return self.secondary_count > 0

    @property
    def has_tertiary_amine(self) -> bool:
        return self.tertiary_count > 0

    @property
    def has_nitrosatable_amine(self) -> bool:
        return self.has_secondary_amine or self.has_tertiary_amine


@dataclass(frozen=True)
class ReagentProfile:
    """시약 화학적 성상 분류."""
    smiles: str
    roles: tuple[str, ...]
    is_nitrite_source: bool
    role_labels_ko: tuple[str, ...]


@dataclass
class ImpurityResult:
    smiles: str
    classification: ImpurityClassification
    rule_id: str
    rule_label: str
    molecular_formula: str
    molecular_weight: float
    db_name: str | None = None
    db_cid: str | None = None
    ashby_tennant_alerts: list[str] = field(default_factory=list)
    highlight_atoms: list[int] = field(default_factory=list)


def _mol_from_smiles(smiles: str) -> Chem.Mol | None:
    if not smiles or not isinstance(smiles, str):
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


def canonical_smiles(smiles: str) -> str | None:
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def _count_substruct(mol: Chem.Mol, smarts: str) -> int:
    try:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            return 0
        return len(mol.GetSubstructMatches(pat))
    except Exception:
        return 0


def detect_secondary_amine_sites(smiles: str) -> AmineSiteInfo:
    return detect_amine_sites(smiles)


def detect_amine_sites(smiles: str) -> AmineSiteInfo:
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return AmineSiteInfo(0, 0)
    sec = _count_substruct(mol, SECONDARY_AMINE_SMARTS)
    ter = _count_substruct(mol, TERTIARY_AMINE_SMARTS)
    return AmineSiteInfo(sec, ter)


def is_nitrite_source_smiles(smiles: str) -> bool:
    """시약 SMILES/텍스트에 아질산염·NOx 패턴이 있는지."""
    if not smiles or not isinstance(smiles, str):
        return False
    low = smiles.lower()
    if any(tok in low for tok in ("nitrite", "nitrous", "nano2", "아질산", "no2", "hon")):
        return True
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return False
    for smarts in NITRITE_SOURCE_SMARTS:
        try:
            pat = Chem.MolFromSmarts(smarts)
            if pat and mol.HasSubstructMatch(pat):
                return True
        except Exception:
            continue
    return False


def classify_reagent(smiles: str) -> ReagentProfile:
    """시약의 산화/환원/알킬화/아질산 성상 분류."""
    roles: list[str] = []
    labels: list[str] = []
    role_label_map = {
        "nitrite": "아질산염 소스 (NDSRI 위험)",
        "oxidant": "산화제",
        "reductant": "환원제",
        "alkylating": "알킬화 시약",
        "acid": "산성 시약",
    }
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return ReagentProfile(smiles=smiles or "", roles=(), is_nitrite_source=False, role_labels_ko=())

    for role, patterns in REAGENT_ROLE_PATTERNS.items():
        for smarts in patterns:
            try:
                pat = Chem.MolFromSmarts(smarts)
                if pat and mol.HasSubstructMatch(pat):
                    if role not in roles:
                        roles.append(role)
                        labels.append(role_label_map.get(role, role))
                    break
            except Exception:
                continue

    nitrite = is_nitrite_source_smiles(smiles) or "nitrite" in roles
    if nitrite and "nitrite" not in roles:
        roles.insert(0, "nitrite")
        labels.insert(0, role_label_map["nitrite"])

    return ReagentProfile(
        smiles=smiles,
        roles=tuple(roles),
        is_nitrite_source=nitrite,
        role_labels_ko=tuple(labels),
    )


def extract_nitrite_mol(reagent_mol: Chem.Mol) -> Chem.Mol | None:
    """시약 분자에서 아질산 이온 부분 구조 추출 (없으면 표준 nitrite 반환)."""
    try:
        for smarts in NITRITE_SOURCE_SMARTS:
            pat = Chem.MolFromSmarts(smarts)
            if pat and reagent_mol.HasSubstructMatch(pat):
                return Chem.MolFromSmiles(NITRITE_SMILES)
    except Exception:
        pass
    return Chem.MolFromSmiles(NITRITE_SMILES)


def _build_reaction(smirks: str) -> rdChemReactions.ChemicalReaction | None:
    try:
        return rdChemReactions.ReactionFromSmarts(smirks)
    except Exception:
        return None


def _sanitize_product_mol(prod: Chem.Mol) -> Chem.Mol | None:
    try:
        if "[H:" in Chem.MolToSmiles(prod):
            prod = Chem.RemoveHs(prod)
        Chem.AssignStereochemistry(prod, force=True, cleanIt=True)
        Chem.SanitizeMol(prod)
        return prod
    except Exception:
        # 2026-09-03 수정: 기존에는 SanitizeMol(catchErrors=True)로 예외만 삼키고 망가진 분자를
        #   그대로 반환했습니다. 예컨대 피롤 N-산화 산물은 kekulize 실패 후 방향족 고리가 조용히
        #   포화 고리로 바뀐 [O-][NH+]1CCCC1(C4H9NO, MW 87.12)이 되어, 날조 구조에 분자식·MW·
        #   M7 등급까지 달렸습니다. sanitize 실패 분자는 불순물로 보고하지 않습니다.
        return None


def _run_smirks_two_component(
    mol_a: Chem.Mol,
    mol_b: Chem.Mol,
    smirks: str,
    *,
    add_h_to_a: bool = False,
) -> list[Chem.Mol]:
    try:
        rxn = _build_reaction(smirks)
        if rxn is None or mol_a is None or mol_b is None:
            return []
        reactant_a = Chem.AddHs(mol_a) if add_h_to_a or "[H:" in smirks else mol_a
        products: list[Chem.Mol] = []

        # 2026-09-03 수정: 반응물 템플릿이 1개인 SMIRKS(방향족 N-산화·카르보닐 환원·
        #   에스터 가수분해)를 항상 2-튜플로 RunReactants에 넘겨 ChemicalParserException이
        #   나고 조용히 0건이 되었습니다. 템플릿 수에 따라 1-튜플/2-튜플을 골라 호출합니다.
        #   1-템플릿일 때 기질(mol_a)에만 적용합니다 — 시약 자신은 기질이 아닙니다.
        if rxn.GetNumReactantTemplates() == 1:
            try:
                outcomes = rxn.RunReactants((reactant_a,))
            except Exception:
                return []
            for outcome in outcomes:
                for prod in outcome:
                    cleaned = _sanitize_product_mol(prod)
                    if cleaned is not None:
                        products.append(cleaned)
            return products

        for order in ((reactant_a, mol_b), (mol_b, reactant_a)):
            try:
                outcomes = rxn.RunReactants(order)
            except Exception:
                continue
            for outcome in outcomes:
                for prod in outcome:
                    cleaned = _sanitize_product_mol(prod)
                    if cleaned is not None:
                        products.append(cleaned)
        return products
    except Exception:
        return []


def _run_smirks_reaction(mol: Chem.Mol, smirks: str) -> list[Chem.Mol]:
    try:
        nitrite = Chem.MolFromSmiles(NITRITE_SMILES)
        if nitrite is None:
            return []
        return _run_smirks_two_component(mol, nitrite, smirks, add_h_to_a="[H:" in smirks)
    except Exception:
        return []


def predict_nitrosation_with_reagent(
    substrate_smiles: str,
    reagent_smiles: str,
) -> tuple[AmineSiteInfo, list[Chem.Mol], list[tuple[str, str, str]]]:
    """
    출발물질(2급/3급 아민) + 그려진 시약(아질산염)으로 NDSRI 예측.
    """
    amine_info = detect_amine_sites(substrate_smiles)
    if not amine_info.has_nitrosatable_amine:
        return amine_info, [], []

    sub_mol = _mol_from_smiles(substrate_smiles)
    reg_mol = _mol_from_smiles(reagent_smiles)
    if sub_mol is None:
        return amine_info, [], []

    nitrite_mol = extract_nitrite_mol(reg_mol) if reg_mol else Chem.MolFromSmiles(NITRITE_SMILES)
    if nitrite_mol is None:
        return amine_info, [], []

    mols: list[Chem.Mol] = []
    meta: list[tuple[str, str, str]] = []
    seen: set[str] = set()

    rules: list[tuple[str, str, str]] = []
    if amine_info.has_secondary_amine:
        rules.append(
            ("reagent_ndsri_secondary", "2급 아민 + 시약 → N-니트로소amine", SECONDARY_AMINE_NITROSATION_SMIRKS)
        )
    if amine_info.has_tertiary_amine:
        rules.append(
            ("reagent_ndsri_tertiary", "3급 아민 + 시약 → 니트로소 탈알킬화 → NDSRI", TERTIARY_AMINE_NITROSATION_SMIRKS)
        )

    for rule_id, rule_label, smirks in rules:
        for prod in _run_smirks_two_component(sub_mol, nitrite_mol, smirks, add_h_to_a="[H:" in smirks):
            try:
                smi = Chem.MolToSmiles(prod, canonical=True)
            except Exception:
                continue
            if not smi or smi in seen:
                continue
            seen.add(smi)
            mols.append(prod)
            meta.append((rule_id, rule_label, smirks))

    return amine_info, mols, meta


def predict_reagent_side_products(
    substrate_smiles: str,
    reagent_smiles: str,
    reagent_profile: ReagentProfile | None = None,
) -> list[tuple[Chem.Mol, str, str]]:
    """
    시약 역할(산화/환원/알킬화/산성)에 따른 부생성물 예측.
    Returns: (mol, rule_id, rule_label)
    """
    sub_mol = _mol_from_smiles(substrate_smiles)
    reg_mol = _mol_from_smiles(reagent_smiles)
    if sub_mol is None or reg_mol is None:
        return []

    profile = reagent_profile or classify_reagent(reagent_smiles)
    if not profile.roles:
        return []

    out: list[tuple[Chem.Mol, str, str]] = []
    seen: set[str] = set()

    for rule_id, rule_label, required_role, smirks in REAGENT_SIDE_REACTIONS:
        if required_role == "nitrite":
            continue
        if required_role not in profile.roles:
            continue
        for prod in _run_smirks_two_component(sub_mol, reg_mol, smirks):
            try:
                smi = Chem.MolToSmiles(prod, canonical=True)
            except Exception:
                continue
            if not smi or smi in seen:
                continue
            seen.add(smi)
            out.append((prod, rule_id, rule_label))

    return out


def predict_reaction_impurities(
    starting_smiles: str,
    *,
    nitrite_present: bool = True,
) -> tuple[AmineSiteInfo, list[Chem.Mol], list[tuple[str, str, str]]]:
    """
    출발물질 SMILES에 SMIRKS 규칙을 적용해 생성 가능한 분자 목록을 반환.

    Returns: (아민정보, 생성 Mol 리스트, (rule_id, rule_label, smirks) 메타)
    """
    amine_info = detect_secondary_amine_sites(starting_smiles)
    if not nitrite_present or not amine_info.has_secondary_amine:
        return amine_info, [], []

    mol = _mol_from_smiles(starting_smiles)
    if mol is None:
        return amine_info, [], []

    mols: list[Chem.Mol] = []
    meta: list[tuple[str, str, str]] = []
    seen: set[str] = set()

    for rule_id, rule_label, smirks in REACTION_REGISTRY:
        for prod in _run_smirks_reaction(mol, smirks):
            smi = Chem.MolToSmiles(prod, canonical=True)
            if smi in seen:
                continue
            seen.add(smi)
            mols.append(prod)
            meta.append((rule_id, rule_label, smirks))

    return amine_info, mols, meta


def build_db_smiles_index(df) -> dict[str, dict]:
    """canonical SMILES → DB 행 메타 (Name, CID)."""
    index: dict[str, dict] = {}
    if df is None or "SMILES" not in df.columns:
        return index
    for _, row in df.iterrows():
        smi = row.get("SMILES")
        if not isinstance(smi, str) or not smi.strip():
            continue
        canon = canonical_smiles(smi)
        if not canon or canon in index:
            continue
        index[canon] = {
            "name": str(row.get("Name", "")),
            "cid": row.get("Compound_CID"),
            "smiles": smi,
        }
    return index


def classify_and_enrich_impurities(
    starting_smiles: str,
    *,
    nitrite_present: bool,
    db_df,
    alert_fn,
) -> tuple[AmineSiteInfo, list[ImpurityResult]]:
    """
    반응 예측 → DB 대조 → Ashby-Tennant 경고 → ImpurityResult 목록.
    """
    amine_info, mols, meta = predict_reaction_impurities(
        starting_smiles, nitrite_present=nitrite_present
    )
    db_index = build_db_smiles_index(db_df)
    results: list[ImpurityResult] = []

    for mol, (rule_id, rule_label, _smirks) in zip(mols, meta):
        smi = Chem.MolToSmiles(mol, canonical=True)
        in_db = smi in db_index
        db_rec = db_index.get(smi, {})
        alerts, highlight_atoms = alert_fn(mol)

        results.append(
            ImpurityResult(
                smiles=smi,
                classification=(
                    ImpurityClassification.VALIDATED
                    if in_db
                    else ImpurityClassification.THEORETICAL
                ),
                rule_id=rule_id,
                rule_label=rule_label,
                molecular_formula=rdMolDescriptors.CalcMolFormula(mol),
                molecular_weight=round(Descriptors.MolWt(mol), 2),
                db_name=db_rec.get("name") if in_db else None,
                db_cid=str(db_rec["cid"]) if in_db and db_rec.get("cid") is not None else None,
                ashby_tennant_alerts=alerts,
                highlight_atoms=highlight_atoms,
            )
        )

    return amine_info, results


def extract_parent_smiles(df, drug_name: str | None = None) -> str | None:
    if df is None or df.empty or "SMILES" not in df.columns:
        return None
    if drug_name and "Name" in df.columns:
        name_lower = drug_name.lower()
        for _, row in df.iterrows():
            val = row.get("Name")
            if isinstance(val, str) and name_lower in val.lower():
                smi = row.get("SMILES")
                if isinstance(smi, str) and smi.strip():
                    return smi.strip()
    first = df["SMILES"].dropna()
    if len(first):
        return str(first.iloc[0]).strip()
    return None


def build_validated_smiles_set(df) -> set[str]:
    return set(build_db_smiles_index(df).keys())
