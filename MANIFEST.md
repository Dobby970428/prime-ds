# 공개본 구성

`_prepare_release.py` 가 155개 파일 (2.4 MB) 를 복사해 만든 폴더입니다.

| 묶음 | 개수 |
|---|---:|
| 엔진 모듈 (p*.py) | 55 |
| 분석 스크립트 (s*.py) | 26 |
| engines\ | 15 |
| webapp\ | 37 |
| templates\ | 1 |
| 배포\ | 7 |
| 루트 부속 파일 | 6 |
| 저장소 문서 | 7 |

## 들어 있지 않은 것

| | 사유 |
|---|---|
| `p45_ionization_probability_filter.py` | 앱도 스크립트도 import 하지 않음. 원고 §2.4 는 물성값이 여과에 쓰이지 않는다고 명시 |
| `p97_build_attentivefp_model.py` | 실제 import 없음(주석 참조만). deepchem 필요, 체크포인트 미배포 |
| `data\CID*.csv` · `drug_index.db` | Search 모듈 전용. 판정 계층은 읽지 않음 |
| `models\` · `cache\` | Ames 데이터셋에서 파생. `data\README.md` 참조 |
| `toxtree_sources\` | Toxtree(GPL-2.0) 배포본 사본 |
| `논문자료\` | 작업 산출물. 상용 도구 출력이 섞여 있어 통째로 제외 |
| `RxnMetrics\` · `_archive\` · `_backup*` · `_selftest_*` | 옛 앱·백업 |

## 실행

```
pip install -r requirements.txt
python p00_flask_app.py
```

`models\ames_ecfp_rf.pkl` 과 `cache\m7_empirical.db` 는 이 폴더에 없습니다.
작업 폴더에서 복사해 넣으면 Layer A·C 가 동작합니다. 두 파일이 없어도 앱은
뜨고, 해당 계층만 비활성으로 경고합니다.
