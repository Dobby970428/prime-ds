# PRIME-DS

**P**redictive **R**isk & **I**mpurity **M**odeling **E**ngine for **D**rug **S**ubstance — an open-source platform that *generates* plausible degradation and process-related impurities from a synthetic route or a formulation and then triages each generated structure under ICH M7(R2).

Most (Q)SAR tools are reactive: they score a structure you give them. PRIME-DS proposes the structures first, then scores them.

> This repository accompanies the article:
>
> Ju, H.; Gong, S.-J.; Yu, S.-H.; Kim, J.-E. *PRIME-DS: An Integrated In Silico Platform for Predicting Degradation and Process-Related Impurities in Drug Substance and Drug Product Development and Supporting ICH M7 Mutagenicity Risk Assessment.* **Pharmaceutics** (submitted).
>
> Archived release: `doi:10.5281/zenodo.23492616`

---

## Quick start

### Option 1 — you already have Python

Python 3.11–3.13. The pinned `scikit-learn` requires ≥ 3.11.

```bash
git clone https://github.com/kookmin-ipl/prime-ds.git
cd prime-ds
pip install -r requirements.txt
python p00_flask_app.py
```

A browser opens at `http://127.0.0.1:8501/`. The server is local and binds to loopback only.

### What leaves the machine

The assessment itself is entirely local — no structure is sent anywhere to be scored. Four things in the interface do reach the network, and all four are optional or replaceable:

| | What | Why it matters |
|---|---|---|
| PubChem REST | the cross-match and the name/structure lookup query `pubchem.ncbi.nlm.nih.gov` | the queried structure or name is sent to NCBI. The assessment layers work without it |
| Ketcher editor | the drawing canvas is an iframe from a **third-party public deployment** of Ketcher | structures drawn in that canvas pass through a host neither we nor Ketcher's authors control. For regulated work, host Ketcher yourself and change the URL in `webapp/views/context.py` and `webapp/static/js/prime.js` |
| cdnjs | the degradation-tree view loads Cytoscape from `cdnjs.cloudflare.com` | no data is sent; the view is inert without it |
| AMBIT bridge | `p52`/`p55` talk to `localhost:8080` if you run the optional AMBIT service | local only |

To run fully offline, skip the PubChem cross-match, self-host Ketcher, and vendor the Cytoscape bundle.

### Option 2 — no Python on the machine

`배포/배포_빌드.bat` (Windows) builds a self-contained folder that carries its own Python runtime. Double-click it once; the result is a portable folder plus a ZIP that runs on any Windows 10+ x64 machine with nothing installed. The build verifies itself before packaging — it loads the model, re-checks the `scikit-learn` version recorded inside the pickle, exercises the alert layer and requests every page. See `배포/README_배포.md`.

### What you need before the assessment layers work

| File | Needed for | If absent |
|---|---|---|
| `models/ames_ecfp_rf.pkl` | Layer C (statistical) | the app starts and warns; Layer C returns nothing |
| `cache/m7_empirical.db` | Layer A (experimental record index) | the app starts and warns; Layer A never hits |
| `drug_index.db` + `data/` | the Search module only | Search opens but returns nothing |

The first two are not in this repository (see **Data availability** below). `requirements.txt` pins `scikit-learn==1.9.0` because the shipped model pickle records that version — a different version either warns and behaves differently or fails to load. Do not relax that pin without retraining.

---

## How the assessment works

Four layers, combined under the **upgrade-only** consensus ICH M7 requires: a statistical prediction may raise a class, never lower one a structural alert has already assigned.

| Layer | What it is | Where | Key parameters |
|---|---|---|---|
| **A** | Experimental record index — Ames and carcinogenicity records looked up by structure | `p60_m7_empirical_db.py`, `cache/m7_empirical.db` | exact structure match; 5,324 records |
| **B** | Expert structural alerts | `p62_toxicity_alerts.py`, `p65_m7_ensemble_matrix.py` | 23 alerts (21 SMARTS + 2 algorithmic, for fused polycyclic aromatics) + 4 metabolic-fragment patterns = 27 patterns; 2 curated exceptions on the primary aromatic amine alert |
| **C** | Statistical (Q)SAR | `p63_deepchem_predictor.py`, `models/ames_ecfp_rf.pkl` | ECFP4, radius 2, 2048 bits; random forest, 100 trees; inference threshold 0.70, escalation threshold 0.85 (`p65`) |
| **D** | Read-across arbitration, fired only when B and C disagree | `p69_read_across_arbitration.py` | Tanimoto ≥ 0.50 against the record index; positive evidence may raise a class, negative evidence never lowers one automatically — it sets the expert-review flag |

Class routing is in `p66_ich_m7.py` (`route_m7_consensus_class_detailed`); the branch order there is what Scheme 1 of the paper draws.

**Applicability domain.** `p73_applicability_domain.py` reports `d_max`, the *maximum* Tanimoto similarity of a query to the training set (higher is closer), with the domain boundary at 0.30. Out-of-domain structures are flagged, not suppressed.

### Generation side

| Module | Contents |
|---|---|
| `p41_degradation_smirks_rules.py` | hydrolysis, oxidation, reduction, photolysis SMIRKS |
| `p30_reaction_engine.py` | secondary- and tertiary-amine nitrosation |
| `engines/p34_name_reaction_process.py`, `engines/p36_suzuki_impurity_rules.py` | named-reaction process impurities, Pd cross-coupling side products |
| `p70`–`p72`, `p18`, `p19` | 143 excipients in nine functional categories, 14 interaction rules |

`s19_rule_inventory.py` counts every rule straight out of the modules and writes the table the paper's Section 2.2 reports — run it rather than counting by hand.

---

## Repository layout

```
p00_flask_app.py          entry point (Flask, 127.0.0.1:8501)
p01 … p96                 engine modules, grouped by number
engines/                  process-impurity and report engines
webapp/                   Flask blueprints, services, templates, i18n (EN/KO)
templates/                read-across report template
s1 … s26                  analysis and figure scripts (see below)
배포/                      portable-build kit (Windows)
data/                     split manifest + how to rebuild the datasets
requirements.txt          development dependencies
```

---

## Reproducing the paper

Every number and figure in the article is the output of a script in this repository. The scripts write into `논문자료/` (paper materials), which is not published here; create it or redirect the paths.

| Script | Produces | Paper |
|---|---|---|
| `s1_setup_layer_c_ames.py` | the scaffold split, the Layer C model | **do not run** — see the warning below |
| `s5_run_pipeline_validation.py` | `pipeline_validation_results.csv`, `report_ames.txt` | Section 3.5, per-compound calls |
| `s11_stats_ad_analysis.py` | `statistics_and_ad_report.txt`, `statistics_and_ad.json`, `ad_tanimoto_validation.csv` | Table 5, Figures 3–4, Section 3.6 |
| `s12_replay_validation_fixed_engine.py` | `replay_fixed_engine_results.csv` | engine-change audit (model held fixed) |
| `s13_export_for_external_tools.py` | blinded `.sdf` / `.smi` + answer key held back | the identical-input premise of Tables 2–4 |
| `s14_score_external_tools.py` | `scaffold_comparison_report.txt`, `scaffold_comparison_merged.csv` | Tables 2 and 4 |
| `s18_ablation_paired_tests.py` | `ablation_paired_tests_report.txt` | Table 3 paired tests |
| `s7_run_case_study.py` | `case_study_report.txt` | Table 7 — NDSRI formation, four APIs |
| `s8_run_excipient_case_study.py` | `excipient_case_study_report.txt` | Section 3.9 — whole-formulation screening |
| `s9_run_binary_compatibility_study.py` | `binary_compatibility_pairs.csv`, `binary_vs_matrix_comparison.csv` | Table 8 — 1:1 binary vs. full matrix |
| `s21_run_ds_process_case_study.py` | `ds_process_case_study_report.txt` | Table 6 — two steps of a losartan route |
| `s22_generation_benefit.py` | `generation_benefit_report.txt` | Tables 9–10, Section 3.12 |
| `s19_rule_inventory.py` | `rule_inventory.csv` / `.txt` | the rule counts in Section 2.2 |
| `s23_make_figures.py` | ROC / applicability-domain panels | Figures 3 and 4 |
| `s24_make_schemes.py` | hydrolysis, nitrosation, oxidation schemes | Schemes 2–4 |
| `s25_make_diagrams.py` | architecture diagram, ICH M7 decision tree | Figure 1, Scheme 1 |
| `s26_build_supplementary.py` | `PRIME-DS_Supplementary_Tables.xlsx` | Tables S1–S4 |

`s23` writes files named `Figure2_*` and `Figure3_*`; these became **Figures 3 and 4** in the final manuscript, after the interface screenshot took the Figure 2 slot.

Schemes 2–4 are drawn by executing the SMIRKS strings that the engine actually holds, and Scheme 1 is traced from the live branch order in `p66`. The figures therefore cannot drift from the code.

> ### Do not run `s1_setup_layer_c_ames.py`
>
> It regenerates the Murcko scaffold split and retrains Layer C. Every comparison in the paper — the validation set, the model, the external-tool scores, the alert exceptions selected on the training split — is tied to the split produced on 2026-09-14. Re-running it silently invalidates all of them. Use the released model and the split manifest in `data/` instead.

### Validation set

Scaffold-separated, 1,288 compounds (588 Ames-positive, 700 negative), no compound sharing an identical ECFP4 fingerprint with any of the 5,309 training compounds.

| | Sensitivity | Specificity | Concordance | MCC |
|---|---|---|---|---|
| Layer B alone (expert rules) | 74.5% | 69.0% | 71.5% | 0.433 |
| Layer C alone (statistical) | 71.1% | 81.0% | 76.5% | 0.524 |
| Four-layer ensemble | 80.8% | 67.0% | 73.3% | 0.478 |

Ensemble specificity is *below* rule-layer specificity by construction: upgrade-only consensus passes every rule-layer false positive through unchanged. That is a property of the ICH M7 combination rule, not a defect of the implementation.

---

## Data availability

**In this repository**

- `data/ames_split_manifest.csv` — 6,597 rows: InChIKey, first-14 InChIKey (connectivity layer), label, and `train` / `validation` assignment. Join it to the source dataset to rebuild the exact splits. See `data/README.md`.

**Not in this repository, and why**

| | Why |
|---|---|
| Ames structures and labels | redistributed under no stated licence by the upstream compilation; the manifest above lets you rebuild the splits from the source |
| `models/ames_ecfp_rf.pkl`, `cache/m7_empirical.db` | derived from that dataset |
| Per-compound Derek Nexus and Sarah Nexus predictions | generated under a Lhasa Limited licence. The paper reports them in aggregate only; per-compound outputs are available from the corresponding author on reasonable request |
| `data/CID*.csv`, `drug_index.db` | a 16,823-file PubChem record cache backing the Search module only. Not used by any assessment layer |
| `toxtree_sources/` | an unpacked copy of the Toxtree distribution (GPL-2.0). Not imported by any code here |

---

## Third-party material and attribution

**Toxtree (GPL-2.0).** Three SMARTS patterns in `p62_toxicity_alerts.py` — the azide/triazene alert and the acyl halide alert — were transcribed verbatim from the Toxtree 3.1.0 source (`mutant/rules/`). Each carries a source comment at its definition. Three further patterns (quinone para/ortho, epoxide/aziridine, α,β-unsaturated lactone) are narrowed variants of Toxtree patterns, and the two fused-polycyclic-aromatic alerts are an independent RDKit reimplementation of the Toxtree decision logic. Toxtree is © Ideaconsult Ltd. and contributors. The unpacked Toxtree distribution is deliberately excluded from this repository; running Toxtree as a benchmark, as the paper does, is use rather than distribution.

**Ames dataset.** Compiled by Therapeutics Data Commons (AMES, single-prediction toxicity), hosted on Harvard Dataverse, `doi:10.7910/DVN/21LKWG`. TDC's primary citation is Xu, C. et al. *J. Chem. Inf. Model.* **2012**. The TDC toolkit is MIT-licensed; the dataset itself carries no stated licence, which is why only the manifest is redistributed here.

**Record index.** 5,274 records from the TDC Ames training split, plus 29 from the Hansen benchmark (Hansen, K. et al. *J. Chem. Inf. Model.* **2009**, *49*, 2077–2081), 12 from Tox21, 8 from ISSCAN.

**Excipient library.** Names follow USP–NF and Ph. Eur. monograph headings; structures are PubChem canonical SMILES; 29 entries derive from the FDA Inactive Ingredient Database (public domain); incompatibility annotations summarise the corresponding entries of the *Handbook of Pharmaceutical Excipients* (Pharmaceutical Press) and are summaries of fact, not reproduced text.

**Libraries.** RDKit (BSD-3-Clause), scikit-learn (BSD-3-Clause), NumPy, pandas, Flask, Pillow, openpyxl, python-docx, requests, rdchiral (all permissive); fpdf2 (LGPL-3.0) and CairoSVG (LGPL-3.0-or-later) are used unmodified as libraries. The structure editor is Ketcher (Apache-2.0, © EPAM Systems), embedded from a third-party public deployment rather than bundled — see "What leaves the machine" above. Cytoscape.js (MIT) is loaded from cdnjs for the degradation-tree view.

**Derek Nexus and Sarah Nexus** are products of Lhasa Limited; **Toxtree** of Ideaconsult Ltd. Neither organisation endorses this work.

---

## Licence

MIT — see [`LICENSE`](LICENSE). The attribution above is part of the licence conditions of the material it names.

## Citing

See [`CITATION.cff`](CITATION.cff), or use the article citation at the top of this file. Cite the archived Zenodo release rather than a branch, so that the version you used stays retrievable.

## Contact

Joo-Eun Kim (corresponding author) — see the article for contact details.
Department of Biopharmaceutical Chemistry / Pharmaceutical Engineering, Kookmin University, Seoul 02707, Republic of Korea.
