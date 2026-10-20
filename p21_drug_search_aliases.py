"""한글 상품명·일반명 → DB 영문 검색어 확장."""

# 한글 -> 영문 변환 (DB가 영문)

from __future__ import annotations

import re

from p22_drug_data import normalize_name

# 한글(또는 별칭) → 영문 INN / 일반명 검색어 목록
_ALIASES: dict[str, list[str]] = {
    "메트포민": ["metformin"],
    "글루코파지": ["metformin"],
    "아스피린": ["aspirin"],
    "아세트아미노펜": ["acetaminophen", "paracetamol"],
    # 2026-09-03 수정: 중복 dict 키 "타이레놀" 한 줄 삭제(값이 동일해 동작은 무해).
    "타이레놀": ["acetaminophen", "paracetamol"],
    "이부프로펜": ["ibuprofen"],
    "부루펜": ["ibuprofen"],
    "아목시실린": ["amoxicillin"],
    "세프트리악손": ["ceftriaxone"],
    "로수바스타틴": ["rosuvastatin"],
    "아토르바스타틴": ["atorvastatin"],
    "리피토": ["atorvastatin"],
    "올메사르탄": ["olmesartan"],
    "암로디핀": ["amlodipine"],
    "로사르탄": ["losartan"],
    "클로피도그렐": ["clopidogrel"],
    "플라빅스": ["clopidogrel"],
    "오메프라졸": ["omeprazole"],
    "란소프라졸": ["lansoprazole"],
    "판토프라졸": ["pantoprazole"],
    "레보티록신": ["levothyroxine"],
    "메트포르민": ["metformin"],
    "시메티딘": ["cimetidine"],
    "디곡신": ["digoxin"],
    "푸로세미드": ["furosemide"],
    "라식스": ["furosemide"],
    "스피로노락톤": ["spironolactone"],
    "와파린": ["warfarin"],
    "헤파린": ["heparin"],
    "인슐린": ["insulin"],
    "글리메피리드": ["glimepiride"],
    "글리벤클라미드": ["glibenclamide"],
    "피오글리타존": ["pioglitazone"],
    "시타글립틴": ["sitagliptin"],
    "엠파글리플로진": ["empagliflozin"],
    "가바펜틴": ["gabapentin"],
    "프레가발린": ["pregabalin"],
    "트라마돌": ["tramadol"],
    "모르핀": ["morphine"],
    "코데인": ["codeine"],
    "디아제팜": ["diazepam"],
    "알프라졸람": ["alprazolam"],
    "졸피뎀": ["zolpidem"],
    "세티리진": ["cetirizine"],
    "로라타딘": ["loratadine"],
    "덱스아메타손": ["dexamethasone"],
    "프레드니솔론": ["prednisolone"],
    "실데나필": ["sildenafil"],
    "타다라필": ["tadalafil"],
    "테오필린": ["theophylline"],
    "살부타몰": ["salbutamol", "albuterol"],
    "벤조일퍼옥사이드": ["benzoyl peroxide"],
    "케토코나졸": ["ketoconazole"],
    "플루코나졸": ["fluconazole"],
    "시프로플록사신": ["ciprofloxacin"],
    "레보플록사신": ["levofloxacin"],
    "아지스로마이신": ["azithromycin"],
    "독시사이클린": ["doxycycline"],
    "메트트렉세이트": ["methotrexate"],
    "리튬": ["lithium"],
    "발프로산": ["valproic acid", "valproate"],
    "카르바마제핀": ["carbamazepine"],
    "페니토인": ["phenytoin"],
    "레보도파": ["levodopa"],
    "도네페질": ["donepezil"],
    "리세드론산": ["risedronate"],
    "알렌드론산": ["alendronate"],
    "카페인": ["caffeine"],
    "니코틴": ["nicotine"],
    "에탄올": ["ethanol"],
}


def has_hangul(text: str) -> bool:
    return any("\uac00" <= c <= "\ud7a3" for c in text)


def expand_search_queries(query: str) -> list[str]:
    """원문 + 별칭으로 확장된 검색어 목록 (중복 제거, 원문 우선)."""
    q = query.strip()
    if not q:
        return []

    out: list[str] = [q]
    qn = normalize_name(q)

    for alias, targets in _ALIASES.items():
        an = normalize_name(alias)
        if not an:
            continue
        if qn == an or an in qn or qn in an:
            out.extend(targets)

    # CID14219 / cid 14219 형식 정규화
    cid_m = re.match(r"^cid\s*(\d+)$", q, re.I)
    if cid_m:
        out.append(f"CID{cid_m.group(1)}")

    seen: set[str] = set()
    unique: list[str] = []
    for item in out:
        key = item.lower()
        if key not in seen:
            seen.add(key)
            unique.append(item)
    return unique
