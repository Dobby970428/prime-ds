# Datasets

This directory holds what can be redistributed. Structures and labels cannot be, so they are not here; the manifest below lets you rebuild them exactly.

## `ames_split_manifest.csv`

The train/validation assignment used throughout the article. 6,597 rows.

| Column | Meaning |
|---|---|
| `inchikey14` | first 14 characters of the InChIKey — the connectivity layer, which is the key the de-duplication step used |
| `inchikey` | full standard InChIKey |
| `label` | Ames outcome: `1` mutagenic, `0` non-mutagenic |
| `split` | `train` or `validation` |

| Split | Compounds | Positive | Negative |
|---|---:|---:|---:|
| `train` | 5,309 | 2,999 | 2,310 |
| `validation` | 1,288 | 588 | 700 |
| total | 6,597 | 3,587 | 3,010 |

Every `inchikey14` is unique across the file, so the two splits share no compound at the connectivity level. The partition itself is by Murcko scaffold, so they also share no scaffold.

## Rebuilding the structures

1. Obtain the Ames dataset from Therapeutics Data Commons — single-prediction toxicity, `AMES`. It is hosted on Harvard Dataverse, `doi:10.7910/DVN/21LKWG`, and TDC's own loader will fetch it:

   ```python
   from tdc.single_pred import Tox
   data = Tox(name='AMES').get_data()      # columns: Drug_ID, Drug, Y
   ```

2. Standardise each SMILES the way the platform does — salt stripping, charge neutralisation, tautomer canonicalisation — with `p74_structure_standardize.py`, then compute the standard InChIKey with RDKit.

3. Join on `inchikey14`. Rows the manifest does not list were dropped by de-duplication, by the removal of scaffolds carrying conflicting labels, or by the exclusions described next.

```python
from rdkit import Chem
import pandas as pd

man = pd.read_csv('data/ames_split_manifest.csv')
# df: your standardised frame with a 'smiles' column
df['inchikey14'] = [
    Chem.MolToInchiKey(Chem.MolFromSmiles(s))[:14] for s in df.smiles
]
merged = df.merge(man, on='inchikey14', how='inner')
assert len(merged) == 6597, len(merged)
```

## How the splits were derived

Starting from 7,278 records in the source compilation:

| Step | Remaining |
|---|---:|
| standardisation, de-duplication on `inchikey14`, removal of scaffolds with conflicting labels | 6,636 |
| 80:20 Murcko scaffold split, no scaffold in both partitions | 5,309 train / 1,327 evaluation |
| drop 27 evaluation compounds whose ECFP4 fingerprint was identical to a training compound | 1,300 |
| drop 12 already held in the experimental record index (Layer A would answer from a record, not a prediction) | **1,288** |

The last two steps matter for interpreting the benchmark: no validation compound shares an identical fingerprint with any training compound, and none is answerable from Layer A.

`s1_setup_layer_c_ames.py` is the script that performed all of this. **Do not re-run it** — it draws a fresh split and retrains the model, which detaches every number in the article from the code that produced it. The manifest exists precisely so that you never have to.

## Not here

| | Why |
|---|---|
| Ames structures and labels as files | the upstream compilation states no licence |
| `models/ames_ecfp_rf.pkl` (Layer C model) | derived from that dataset |
| `cache/m7_empirical.db` (Layer A record index, 5,324 rows) | 5,274 of its rows are the TDC training split |
| `data/CID*.csv`, `drug_index.db` | a PubChem record cache used only by the Search module; no assessment layer reads it |
| Per-compound Derek Nexus / Sarah Nexus predictions | generated under a Lhasa Limited licence; the article reports them in aggregate, and per-compound outputs are available from the corresponding author on reasonable request |
| Per-compound Toxtree predictions | may be released on request; omitted here only for size |

## Citing the dataset

Cite both the compilation and its primary source:

- Huang, K. et al. Therapeutics Data Commons. `doi:10.7910/DVN/21LKWG`
- Xu, C. et al. In Silico Prediction of Chemical Ames Mutagenicity. *J. Chem. Inf. Model.* **2012**.

The 29 Hansen-benchmark records in the Layer A index come from Hansen, K. et al. *J. Chem. Inf. Model.* **2009**, *49*, 2077–2081.
