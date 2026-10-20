"""
구조 표준화 — 적용범위(AD) 산출의 단일 기준.

왜 있는가
    s1_setup_layer_c_ames.py 는 2026-09-03 수정(C-1)부터 학습셋·검증셋을
    표준화한 뒤 저장한다 (s1 395-402행: 염 제거 → 전하 중화 → 호변이성체
    정규화 → canonical SMILES).

    그런데 s11 과 p73 의 fingerprint() 는 Chem.MolFromSmiles() 만 부르고
    표준화를 하지 않았다. 학습셋은 s1 이 이미 표준화해 두었으므로 문제가
    드러나지 않았지만, 사례연구 생성 구조(SMIRKS 산물)는 표준화되지 않은
    원본이다. 즉 **표준화된 학습셋과 표준화되지 않은 질의 구조로 Tanimoto
    를 재고 있었다.**

    구 분할(랜덤)에서는 학습셋도 표준화하지 않았으므로 원본 대 원본으로
    일관됐다. 스캐폴드 분할 재실행(2026-09-14)이 이 불일치를 만들었다.

    이 모듈이 표준화 정의를 한 곳에 두고, s11 과 p73 이 함께 쓴다.

정의 — s1 의 standardize_mol 과 **글자 그대로 동일해야 한다**
    rdMolStandardize.LargestFragmentChooser().choose()   염·짝이온 제거
    rdMolStandardize.Uncharger().uncharge()              전하 중화
    rdMolStandardize.TautomerEnumerator().Canonicalize() 호변이성체 정규화
    Chem.SanitizeMol()

    실패하면 원본 구조를 그대로 쓴다 (s1 398행과 같은 처리).

비용
    TautomerEnumerator.Canonicalize 가 느리다. 구조당 수십 밀리초 수준이므로
    수천 건을 처리하면 분 단위가 걸린다. p73 은 학습셋 지문을 캐시하므로
    최초 1회만 부담한다. s11 은 실행할 때마다 부담한다.

단독 실행
    python p74_structure_standardize.py      자기 점검
"""

from __future__ import annotations

from typing import Any

_CACHE: dict[str, Any] = {}

# 표준화 정의가 바뀌면 이 값을 올린다. p73 의 캐시 서명에 들어가므로
# 값이 바뀌면 학습셋 지문 캐시가 자동으로 폐기된다.
STANDARDIZER_VERSION = 1


def _standardizers():
    """rdMolStandardize 도구를 한 번만 만들어 재사용한다. (s1 _standardizers 와 동일)"""
    if "ready" not in _CACHE:
        from rdkit.Chem.MolStandardize import rdMolStandardize

        _CACHE["lfc"] = rdMolStandardize.LargestFragmentChooser()
        _CACHE["unc"] = rdMolStandardize.Uncharger()
        _CACHE["taut"] = rdMolStandardize.TautomerEnumerator()
        _CACHE["ready"] = True
    return _CACHE["lfc"], _CACHE["unc"], _CACHE["taut"]


def standardize_mol(mol):
    """염 제거 → 전하 중화 → 호변이성체 정규화. 실패하면 None. (s1 과 동일)"""
    from rdkit import Chem

    lfc, unc, taut = _standardizers()
    try:
        m = lfc.choose(mol)
        m = unc.uncharge(m)
        m = taut.Canonicalize(m)
        Chem.SanitizeMol(m)
        return m
    except Exception:  # noqa: BLE001
        return None


def standardized_mol_from_smiles(smiles: str):
    """
    SMILES → 표준화된 Mol. 파싱 실패면 None.

    표준화에 실패하면 s1 398행과 같이 **원본 구조를 그대로** 돌려준다.
    지문 산출이 통째로 실패하는 것보다 낫고, s1 이 학습셋에 적용한 처리와
    같아야 하기 때문이다.
    """
    from rdkit import Chem

    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None
    std = standardize_mol(mol)
    return std if std is not None else mol


def _self_check() -> int:
    """표준화가 실제로 작동하는지 몇 가지로 확인한다."""
    from rdkit import Chem

    cases = [
        ("염 제거", "CC(=O)Oc1ccccc1C(=O)[O-].[Na+]", "CC(=O)Oc1ccccc1C(=O)O"),
        ("전하 중화", "CC(=O)[O-]", "CC(=O)O"),
        ("이미 표준", "CN(C)N=O", "CN(C)N=O"),
    ]

    print()
    print("=" * 62)
    print("  p74 구조 표준화 자기 점검")
    print("=" * 62)
    print()

    failed = 0
    for label, smi, expect in cases:
        mol = standardized_mol_from_smiles(smi)
        got = Chem.MolToSmiles(mol) if mol is not None else "(실패)"
        exp = Chem.MolToSmiles(Chem.MolFromSmiles(expect))
        ok = got == exp
        failed += 0 if ok else 1
        print(f"  [{'O' if ok else 'X'}] {label:<10} {smi}")
        print(f"      → {got}")
        if not ok:
            print(f"      기대 {exp}")
    print()

    # 호변이성체 — 두 표기가 같은 구조로 모이는지
    a = standardized_mol_from_smiles("O=C1CCCCC1")
    b = standardized_mol_from_smiles("OC1=CCCCC1")
    if a is not None and b is not None:
        same = Chem.MolToSmiles(a) == Chem.MolToSmiles(b)
        print(f"  [{'O' if same else 'X'}] 호변이성체   케토/엔올이 같은 구조로 정규화")
        print(f"      케토 → {Chem.MolToSmiles(a)}")
        print(f"      엔올 → {Chem.MolToSmiles(b)}")
        failed += 0 if same else 1
    print()

    print(f"  실패 {failed}건")
    print()
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_self_check())
