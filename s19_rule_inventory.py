"""
반응 규칙 전수 조사 — 원고 2.2절 표의 근거.

왜 있는가
    논문 초안 2.2 절의 규칙 수 표(가수분해 6 · 산화 7 · 환원 2 · 광분해 2 ·
    니트로소화 2 · Name Reaction 33 · 부형제 2)는 손으로 적은 것이며
    코드와 맞지 않는다. 예컨대 p41 에 실제로 있는 규칙은 7종이고,
    부형제 규칙은 2종이 아니라 14종이다.

    규칙이 여러 모듈에 흩어져 있어 눈으로 세면 또 틀린다. 이 스크립트가
    모듈을 직접 읽어 표를 만든다. 원고의 표는 여기 출력을 옮겨 적는다.

어떻게 세는가
    각 모듈의 최상위 이름과 **원본 코드**를 훑어 '>>' 가 들어 있는 문자열을 담은
    컨테이너(dict · tuple · list)를 찾는다. SMIRKS 변수명을 미리 알 필요가
    없으므로 모듈이 바뀌어도 따라간다.

    찾은 SMIRKS 는 RDKit 으로 파싱해 보고, 실패한 것과 반응물·생성물
    템플릿 수를 함께 보고한다. 파싱 실패는 «조용히 빠지는» 규칙이므로
    원고에 세면 안 된다.

출력
    논문자료/02_검증결과/rule_inventory.txt
    논문자료/02_검증결과/rule_inventory.csv

사용법
    py s19_rule_inventory.py
"""

from __future__ import annotations

import importlib
from pathlib import Path

APP = Path(__file__).resolve().parent
OUT_DIR = APP / "논문자료" / "02_검증결과"
OUT_TXT = OUT_DIR / "rule_inventory.txt"
OUT_CSV = OUT_DIR / "rule_inventory.csv"

# (모듈 경로, 원고에서 부를 이름)
MODULES = [
    ("p18_excipient_reaction_rules", "부형제 상호작용"),
    ("p30_reaction_engine", "니트로소화 · 시약 부반응"),
    ("p31_process_reaction_rules", "공정 반응 규칙"),
    ("p32_process_reaction_engine", "공정 반응 엔진"),
    ("p41_degradation_smirks_rules", "분해 (SMIRKS 카탈로그)"),
    ("p42_degradation_rag_engine", "분해 (RAG 엔진)"),
    ("engines.p34_name_reaction_process", "Name reaction"),
    ("engines.p35_universal_reaction_impurity_engine", "공정 불순물 (범용)"),
    ("engines.p36_suzuki_impurity_rules", "Suzuki 불순물"),
]

_log: list[str] = []


def w(msg: str = "") -> None:
    print(msg)
    _log.append(msg)


# 규칙 이름이 아니라 «필드 이름» 인 키들. 이게 라벨로 잡히면
# main_smarts · smirks 같은 쓸모없는 이름이 표에 들어간다.
_FIELD_KEYS = frozenset({
    "smirks", "smarts", "main_smirks", "side_smirks", "main_smarts",
    "side_smarts", "rxn", "reaction", "reaction_smarts", "template",
})

# 실행 경로에 없는 컨테이너. p30 의 RETIRED_* 는 «어디에서도 실행되지
# 않습니다» 라고 코드에 적혀 있으므로 세면 안 된다.
_SKIP_PREFIX = ("RETIRED_", "_RETIRED", "DEPRECATED_", "OLD_", "LEGACY_")

# 반응 템플릿이 아니라 **문자열 부분일치용 차단 목록**인 컨테이너.
# p41._METABOLIC_SMIRKS_FRAGMENTS 는 «이 SMIRKS 가 대사 계열인가» 를
# 판정하려고 일부러 잘린 조각을 담고 있다. '>>' 가 들어 있다는 이유로
# 템플릿으로 세면 규칙 수가 부풀고, 파싱 실패가 «고장» 처럼 보인다.
# (2026-09-15 — 실제로 그렇게 오해했다.)
_SKIP_EXACT = frozenset({
    "_METABOLIC_SMIRKS_FRAGMENTS",
    "_METABOLIC_RULE_NAME_BLOCKLIST",
})

# 딕셔너리 안에서 규칙 이름으로 쓸 만한 키.
_NAME_KEYS = ("name", "rule", "rule_name", "id", "label", "title", "reaction_name")


def _is_smirks(v: object) -> bool:
    return isinstance(v, str) and ">>" in v and len(v) < 4000


def _dict_label(d: dict) -> str:
    """딕셔너리가 스스로 들고 있는 이름을 꺼낸다."""
    for k in _NAME_KEYS:
        v = d.get(k)
        if isinstance(v, str) and v.strip():
            return v.strip()
    return ""


def _walk(value: object, depth: int = 0, inherited: str = "") -> list[tuple[str, str]]:
    """
    컨테이너를 훑어 (규칙 이름, SMIRKS) 쌍을 모은다.

    이름은 **바깥에서 안으로** 물려준다. 안쪽 필드 키(smirks 등)가
    이름으로 잡히면 표에 쓸 수 없기 때문이다.
    """
    found: list[tuple[str, str]] = []
    if depth > 4:
        return found

    if _is_smirks(value):
        found.append((inherited, str(value)))

    elif isinstance(value, dict):
        own = _dict_label(value)
        for k, v in value.items():
            key = str(k)
            if key in _NAME_KEYS:
                continue
            # 바깥 키가 규칙 이름이면 그것을 쓰고, 필드 이름이면 물려받은
            # 이름을 유지한다.
            label = inherited or own
            if key not in _FIELD_KEYS and not key.isdigit():
                label = key if not own else own
            found.extend(_walk(v, depth + 1, label))

    elif isinstance(value, (list, tuple, set, frozenset)):
        for i, item in enumerate(value):
            if (isinstance(item, (list, tuple)) and len(item) >= 2
                    and isinstance(item[0], str) and not _is_smirks(item[0])):
                # (이름, SMIRKS, ...) 형태
                for sub in item[1:]:
                    found.extend(_walk(sub, depth + 1, item[0]))
            else:
                found.extend(_walk(item, depth + 1, inherited))

    return found


def _source_scan(mod) -> list[tuple[str, str]]:
    """
    최상위 이름으로 노출되지 않는 템플릿을 **원본 코드**에서 찾는다.

    2026-09-28 추가. 기존 스캔은 `dir(mod)` 로 모듈 최상위 이름만 훑고
    클래스는 건너뛴다(main 의 `callable(value) or isinstance(value, type)`).
    그런데 `engines.p35_universal_reaction_impurity_engine` 는 템플릿을
    클래스 `__init__` 안에서 `self.xxx = AllChem.ReactionFromSmarts("...")`
    로 만든다. 인스턴스 속성이므로 최상위에도 클래스 속성에도 없고,
    그 결과 **7종이 통째로 빠져 전체가 105종으로 집계되었다**(원고 본문은
    112종). 인스턴스를 만들어 보는 대신 AST 로 문자열 상수를 읽는다 —
    부작용이 없고 import 만 되면 동작한다.

    반환: (이름, SMIRKS) 목록. 이름은 대입 대상에서 딴다
          (`self.cross_coupling_rxn` → `cross_coupling_rxn`).
    """
    import ast
    import inspect

    try:
        src = inspect.getsource(mod)
    except (OSError, TypeError):
        return []
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return []

    out: list[tuple[str, str]] = []
    for node in ast.walk(tree):
        # 주석형 대입(`X: dict = {...}`)은 ast.AnnAssign 이라 함께 본다.
        if isinstance(node, ast.Assign) and node.targets:
            tgt = node.targets[0]
        elif isinstance(node, ast.AnnAssign) and node.value is not None:
            tgt = node.target
        else:
            continue
        if isinstance(tgt, ast.Attribute):
            name = tgt.attr
        elif isinstance(tgt, ast.Name):
            name = tgt.id
        else:
            continue
        if name in _SKIP_EXACT or name.startswith(_SKIP_PREFIX):
            continue
        for sub in ast.walk(node.value):
            if isinstance(sub, ast.Constant) and _is_smirks(sub.value):
                out.append((name, str(sub.value)))
    return out


def main() -> int:
    from rdkit import RDLogger
    from rdkit.Chem import AllChem

    RDLogger.DisableLog("rdApp.*")

    w("=" * 78)
    w("  반응 규칙 전수 조사 — 원고 2.2 절 표의 근거")
    w("=" * 78)
    w()
    w("  각 모듈의 최상위 이름에서 '>>' 를 포함한 문자열을 찾고,")
    w("  최상위에 노출되지 않는 것(클래스 내부)은 원본 코드에서 보충한다.")
    w("  RDKit 파싱에 실패한 규칙은 실행 시 조용히 빠지므로 따로 표시한다.")
    w()

    rows: list[dict] = []
    total_ok = total_bad = 0

    for mod_path, label in MODULES:
        try:
            mod = importlib.import_module(mod_path)
        except Exception as exc:  # noqa: BLE001
            w(f"  [건너뜀] {mod_path} — import 실패: {type(exc).__name__}: {exc}")
            continue

        seen: set[str] = set()
        entries: list[tuple[str, str, str]] = []   # (컨테이너, 이름, smirks)

        skipped_retired = 0
        skipped_filter = 0
        skipped_empty = 0
        for attr in dir(mod):
            if attr.startswith("__"):
                continue
            if attr in _SKIP_EXACT:
                skipped_filter += 1
                continue
            if attr.startswith(_SKIP_PREFIX):
                # 은퇴한 규칙 문자열 — 실행 경로가 아니므로 세지 않는다.
                try:
                    if any(True for _ in _walk(getattr(mod, attr))):
                        skipped_retired += 1
                except Exception:  # noqa: BLE001
                    pass
                continue
            try:
                value = getattr(mod, attr)
            except Exception:  # noqa: BLE001
                continue
            if callable(value) or isinstance(value, type):
                continue
            for nm, smirks in _walk(value):
                if smirks in seen:
                    continue
                seen.add(smirks)
                entries.append((attr, nm or attr, smirks))

        # 최상위 이름에 없던 것(클래스 내부 등)을 원본 코드에서 보충한다.
        for nm, smirks in _source_scan(mod):
            if smirks in seen:
                continue
            seen.add(smirks)
            entries.append(("(클래스 내부)", nm, smirks))

        if not entries:
            w(f"  [없음] {mod_path} — SMIRKS 를 찾지 못했다")
            continue

        ok = bad = 0
        w("-" * 78)
        w(f"  {label}   ({mod_path})")
        w("-" * 78)
        for container, name, smirks in sorted(entries, key=lambda e: (e[0], e[1])):
            try:
                rxn = AllChem.ReactionFromSmarts(smirks)
            except Exception:  # noqa: BLE001
                rxn = None
            if rxn is None:
                bad += 1
                status = "파싱 실패"
                n_r = n_p = -1
            else:
                n_r = rxn.GetNumReactantTemplates()
                n_p = rxn.GetNumProductTemplates()
                if n_r == 0 or n_p == 0:
                    # 2026-09-29: 반응물이나 생성물이 0 인 것은 반응 템플릿이
                    # 아니다. 원본 스캔을 도입하면서 p34 의 `reactant_templates`
                    # 같은 «'>>' 를 포함하지만 반응이 아닌» 문자열이 잡혀
                    # 전체가 112 → 113 종으로 부풀었다. RDKit 파싱에는 성공하나
                    # 실행 시 아무 산물도 내지 않으므로 세지 않는다.
                    skipped_empty += 1
                    continue
                ok += 1
                status = "정상"
            rows.append({
                "module": mod_path, "group": label, "container": container,
                "rule": name, "status": status,
                "n_reactant": n_r, "n_product": n_p, "smirks": smirks,
            })
            mark = " " if rxn is not None else "✗"
            w(f"   {mark} {name[:44]:<46} {container[:24]:<26} "
              f"{'' if rxn is None else f'{n_r}→{n_p}'}")

        total_ok += ok
        total_bad += bad
        w(f"      소계 {ok + bad}종 (정상 {ok} · 파싱 실패 {bad})")
        if skipped_retired:
            w(f"      은퇴 컨테이너 {skipped_retired}개는 세지 않았다 "
              f"(RETIRED_* · DEPRECATED_* 등)")
        if skipped_filter:
            w(f"      차단 목록 컨테이너 {skipped_filter}개는 세지 않았다 "
              f"(반응 템플릿이 아니라 문자열 부분일치용)")
        w()

    w("=" * 78)
    w("  요약")
    w("=" * 78)
    w()
    w(f"  전체 {total_ok + total_bad}종 — 정상 {total_ok} · 파싱 실패 {total_bad}")
    w()
    w("  ⚠ 같은 SMIRKS 가 여러 모듈에 복제돼 있으면 중복으로 세어진다.")
    w("    모듈별 소계를 그대로 더하지 말고, 원고에는 «어느 모듈이 실행")
    w("    경로에 있는가» 를 확인한 뒤 해당 모듈만 세어 적을 것.")
    w("    (예 — 논문자료/06_재현스크립트/ 의 사본은 실행 경로가 아니다.)")
    w()
    w("  ⚠ 파싱 실패한 규칙은 실행 시 빠지므로 원고 표에 세면 안 된다.")
    w()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_TXT.write_text("\n".join(_log), encoding="utf-8")
    try:
        import pandas as pd
        pd.DataFrame(rows).to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
        print(f"\n목록 : {OUT_CSV}")
    except Exception as exc:  # noqa: BLE001
        print(f"(CSV 저장 실패: {exc})")
    print(f"기록 : {OUT_TXT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
