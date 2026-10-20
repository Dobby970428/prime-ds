"""부형제 데이터의 영문 표시 이름과 주의 문구.

엔진 데이터(`p70`/`p71`/`p72`)는 건드리지 않았습니다. 화면에 영어로
보여줄 때만 이 표를 씁니다. 한국어 화면에서는 엔진의 label_ko와
alerts를 그대로 씁니다.

영문 이름은 USP-NF / Ph.Eur. 표제어를 따랐고, 널리 쓰이는 상품명은
괄호에 남겼습니다.

2026-08-25 작성.
"""

from __future__ import annotations

# key → (English name, English alert)
EXCIPIENT_EN: dict[str, tuple[str, str]] = {
    # --- Binder / Diluent ---
    "Lactose_Monohydrate": (
        "Lactose monohydrate",
        "🔥 Reducing sugar — Maillard browning and condensation impurities with "
        "primary/secondary amine drug substances",
    ),
    "Lactose_Anhydrous": (
        "Lactose, anhydrous",
        "🔥 Open-chain aldehyde persists even in low-moisture formulations — Maillard risk",
    ),
    "Microcrystalline_Cellulose": (
        "Microcrystalline cellulose (MCC / Avicel)",
        "Moisture sorption accelerates solid-state hydrolysis, particularly of ester "
        "and amide drug substances",
    ),
    "Mannitol": (
        "Mannitol",
        "Alcohol dehydration by-products; residual oxidation catalyst on cross-contamination",
    ),
    "Sorbitol": (
        "Sorbitol",
        "Hygroscopic sugar alcohol — Maillard reaction and trace aldehyde formation",
    ),
    "Dibasic_Calcium_Phosphate": (
        "Dibasic calcium phosphate (Emcompress)",
        "Alkaline diluent — local pH rise promotes base-catalysed hydrolysis",
    ),
    "Calcium_Carbonate": (
        "Calcium carbonate",
        "Releases CO₂ and alkali in the formulation — accelerates degradation of "
        "acid-labile drug substances",
    ),
    "Corn_Starch": (
        "Starch, corn",
        "Residual reducing sugars — Maillard reaction and discolouration with amine "
        "drug substances",
    ),
    "Pregelatinized_Starch": (
        "Starch, pregelatinised",
        "High-moisture binder — promotes Maillard reaction and solid-state hydrolysis",
    ),
    "Sucrose": (
        "Sucrose",
        "🔥 Inverts to reducing sugars on thermal or acid degradation — Maillard "
        "reaction and caramelisation",
    ),
    "Isomalt": (
        "Isomalt",
        "Sugar alcohol mixture — reducing residues persist even at low moisture",
    ),

    # --- Lubricant / Glidant ---
    "Magnesium_Stearate": (
        "Magnesium stearate",
        "🔥 Free fatty acid and Mg²⁺ — catalyses oxidation and hydrolysis; reduces dissolution",
    ),
    "Stearic_Acid": (
        "Stearic acid",
        "Secondary impurities via esterification or salt formation with amine and "
        "alcohol drug substances",
    ),
    "Talc": (
        "Talc",
        "Trace transition metals — may initiate radical autoxidation",
    ),
    "Colloidal_Silicon_Dioxide": (
        "Silica, colloidal anhydrous (Aerosil)",
        "High surface area concentrates moisture and oxidants at the interface",
    ),
    "Sodium_Stearyl_Fumarate": (
        "Sodium stearyl fumarate (Pruv)",
        "Fumarate ester — reversible Michael-type addition with amine drug substances",
    ),

    # --- Disintegrant ---
    "Crospovidone": (
        "Crospovidone (Kollidon CL)",
        "🔥 Residual peroxides — N-oxide and sulfoxide formation",
    ),
    "Croscarmellose_Sodium": (
        "Croscarmellose sodium (Ac-Di-Sol)",
        "Rapid swelling causes local pH shifts and accelerates ambient hydrolysis",
    ),
    "Sodium_Starch_Glycolate": (
        "Sodium starch glycolate (Explotab)",
        "Starch-derived reducing sugars plus carboxylic acid — cross Maillard reaction",
    ),
    "Low_Substituted_HPC": (
        "Hydroxypropylcellulose, low-substituted (L-HPC)",
        "Promotes water uptake — accelerates solid-state hydrolysis",
    ),

    # --- Polymer / Film coating ---
    "HPMC": (
        "Hypromellose (HPMC)",
        "Moisture barrier during film coating — core/coat pH imbalance degradation",
    ),
    "Povidone_K30": (
        "Povidone K30 (PVP K30)",
        "🔥 Residual peroxides — promotes oxidative impurities",
    ),
    "Polyethylene_Glycol_4000": (
        "Macrogol 4000 (PEG 4000)",
        "🔥 Formic acid and formaldehyde — formylation and NDSRI precursor",
    ),
    "Ethyl_Cellulose": (
        "Ethylcellulose",
        "Residual organic solvent plus oxidative degradation acids interact with the "
        "drug substance",
    ),
    "Hydroxypropyl_Cellulose": (
        "Hydroxypropylcellulose (Klucel)",
        "Binder and coating agent — moisture sorption causes solid-state instability",
    ),
    "Methacrylic_Acid_Copolymer": (
        "Methacrylic acid copolymer (Eudragit L/S)",
        "pH-dependent dissolution — base-catalysed hydrolysis in intestinal alkaline media",
    ),
    "Poloxamer_407": (
        "Poloxamer 407 (Pluronic F127)",
        "Thermoresponsive gel — residual autoxidation peroxides",
    ),

    # --- Surfactant ---
    "Polysorbate_80": (
        "Polysorbate 80 (Tween 80)",
        "🔥 Unsaturated fatty acid — peroxide formation and epoxidation degradation",
    ),
    "Polysorbate_20": (
        "Polysorbate 20 (Tween 20)",
        "Ester surfactant — hydrolysis and peroxide formation",
    ),
    "Sodium_Lauryl_Sulfate": (
        "Sodium lauryl sulfate (SLS)",
        "Strong surfactant — increases water accessibility and accelerates hydrolysis",
    ),
    "Benzalkonium_Chloride": (
        "Benzalkonium chloride",
        "Quaternary ammonium — ion-pair impurities with anionic drug substances",
    ),

    # --- Buffer / pH modifier ---
    "Citric_Acid": (
        "Citric acid",
        "Acidic pH — promotes acid hydrolysis of esters, amides and β-lactams",
    ),
    "Sodium_Citrate": (
        "Sodium citrate",
        "Alkaline buffer — activates base-catalysed hydrolysis at pH 6–7",
    ),
    "Tartaric_Acid": (
        "Tartaric acid",
        "Dicarboxylic acid — salt formation and solid-state reaction with amine "
        "drug substances",
    ),
    "Sodium_Bicarbonate": (
        "Sodium bicarbonate",
        "Releases CO₂ — local pH rise affects acid-labile drug substances",
    ),
    "Tromethamine_TRIS": (
        "Trometamol (TRIS)",
        "Primary amine buffer — Schiff base formation with aldehyde drug substances",
    ),

    # --- Sweetener / Flavour ---
    "Aspartame": (
        "Aspartame",
        "🔥 Peptide bond — Asp-Phe diketopiperazine forms under moisture and pH stress",
    ),
    "Sucralose": (
        "Sucralose",
        "Chlorinated sugar — releases HCl on thermal degradation, driving acid hydrolysis",
    ),
    "Saccharin_Sodium": (
        "Saccharin sodium",
        "Organic acid salt — affects pH and forms ion pairs",
    ),
    "Peppermint_Oil": (
        "Peppermint oil (representative: menthol)",
        "Terpene — oxidation and extraction of the drug substance",
    ),

    # --- Preservative / Antioxidant ---
    "Methylparaben": (
        "Methyl parahydroxybenzoate",
        "Ester preservative — hydrolyses to residual p-hydroxybenzoic acid",
    ),
    "Propylparaben": (
        "Propyl parahydroxybenzoate",
        "Phenolic ester — alkaline hydrolysis liberates free phenol",
    ),
    "Benzyl_Alcohol": (
        "Benzyl alcohol",
        "Oxidises to benzaldehyde — Maillard reaction and Schiff base formation",
    ),
    "Butylated_Hydroxytoluene_BHT": (
        "Butylated hydroxytoluene (BHT)",
        "Antioxidant — inhibits oxidation but generates quinone-type degradation products",
    ),
    "Butylated_Hydroxyanisole_BHA": (
        "Butylated hydroxyanisole (BHA)",
        "Antioxidant — quinone by-products",
    ),
    "Ascorbic_Acid": (
        "Ascorbic acid",
        "Reducing agent — metal-catalysed redox with the drug substance; supports "
        "Maillard reaction",
    ),
    "Edetate_Disodium_EDTA": (
        "Disodium edetate (EDTA)",
        "Metal chelator — alters degradation pathways by sequestering catalytic metals",
    ),
    "Sodium_Benzoate": (
        "Sodium benzoate",
        "Preservative — liberates benzoic acid at acidic pH",
    ),
    "Potassium_Sorbate": (
        "Potassium sorbate",
        "Unsaturated acid salt — oxidative degradation",
    ),
    "Sodium_Metabisulfite": (
        "Sodium metabisulfite",
        "🔥 Reducing agent — reacts with amine and aldehyde drug substances",
    ),

    # --- Solvent / Co-solvent ---
    "Propylene_Glycol": (
        "Propylene glycol",
        "Oxidises to lactic acid and acetone — local pH shifts",
    ),
    "Glycerol": (
        "Glycerol",
        "Hygroscopic — raises water activity and supports Maillard reaction",
    ),
    "Ethanol": (
        "Ethanol",
        "Oxidises to acetaldehyde — condensation with amine drug substances",
    ),
    "Isopropyl_Alcohol": (
        "Isopropyl alcohol (IPA)",
        "Solvent — oxidises to acetone",
    ),
    "Benzyl_Benzoate": (
        "Benzyl benzoate",
        "Ester solvent — hydrolysis by-products",
    ),
    "Light_Mineral_Oil": (
        "Light mineral oil",
        "Hydrocarbon — extractables and residual solvent",
    ),

    # --- 오픈소스 추가분 (FDA IID) ---
    "Gelatin": (
        "Gelatin",
        "Protein hydrolysis — Maillard reaction with amine drug substances and "
        "isomerisation risk",
    ),
    "Xanthan_Gum": (
        "Xanthan gum",
        "Polysaccharide — retains moisture and promotes solid-state hydrolysis",
    ),
    "Carboxymethylcellulose_Sodium": (
        "Carmellose sodium (CMC)",
        "Viscosity modifier — raises local water activity",
    ),
    "Alginic_Acid": (
        "Alginic acid",
        "Polysaccharide acid — pH buffering and ion exchange",
    ),
    "Acacia_Gum": (
        "Acacia (gum arabic)",
        "Natural polysaccharide — supports Maillard reaction",
    ),
    "Sodium_Chloride": (
        "Sodium chloride",
        "Ionic strength shift — affects solubility and hydrolysis kinetics",
    ),
    "Dextrose": (
        "Dextrose (D-glucose)",
        "🔥 Reducing sugar — Maillard impurities",
    ),
    "Calcium_Phosphate_Tribasic": (
        "Tribasic calcium phosphate",
        "Alkaline diluent — pH rise",
    ),
    "Titanium_Dioxide": (
        "Titanium dioxide",
        "Photocatalytic surface — supports photolytic degradation in film coatings",
    ),
    "Glyceryl_Monostearate": (
        "Glyceryl monostearate",
        "Fatty acid ester — hydrolysis and oxidation",
    ),
    "Hydrogenated_Castor_Oil": (
        "Castor oil, hydrogenated",
        "Lipid — oxidation and peroxide formation",
    ),
    "Hypromellose_Phthalate": (
        "Hypromellose phthalate (HPMCP)",
        "Enteric coating — pH-dependent dissolution",
    ),
    "Cellulose_Acetate_Phthalate": (
        "Cellacefate (CAP)",
        "Enteric coating — residual acetic acid from hydrolysis",
    ),
    "Povidone_K90": (
        "Povidone K90",
        "Polymeric binder — residual peroxides",
    ),
    "Polysorbate_60": (
        "Polysorbate 60 (Tween 60)",
        "Ester surfactant — peroxide formation",
    ),
    "Lecithin": (
        "Lecithin",
        "Phospholipid — oxidation and radical formation",
    ),
    "Sodium_Phosphate_Dibasic": (
        "Disodium phosphate",
        "Buffer salt — base-catalysed hydrolysis at pH 7–8",
    ),
    "Potassium_Phosphate": (
        "Potassium phosphate",
        "Buffer — ionic strength and pH shifts",
    ),
    "Sodium_Hydroxide": (
        "Sodium hydroxide",
        "Strong base — promotes base-catalysed hydrolysis",
    ),
    "Hydrochloric_Acid": (
        "Hydrochloric acid",
        "Strong acid — promotes acid hydrolysis",
    ),

    # ======================================================================
    # 2026-08-25 확장분 (p19_excipient_library_ext)
    # ======================================================================

    # --- Binder / Diluent ---
    "Maltodextrin": ("Maltodextrin",
        "🔥 Reducing end group — Maillard reaction with primary amines"),
    "Maltose": ("Maltose", "🔥 Reducing sugar — Maillard browning"),
    "Fructose": ("Fructose", "🔥 Ketose — Maillard reaction and carbonyl chemistry"),
    "Trehalose": ("Trehalose",
        "Non-reducing sugar — low Maillard risk; lyophilisation stabiliser"),
    "Xylitol": ("Xylitol", "Sugar alcohol — non-reducing; hygroscopic"),
    "Erythritol": ("Erythritol", "Sugar alcohol — non-reducing"),
    "Calcium_Sulfate_Dihydrate": ("Calcium sulfate dihydrate",
        "Releases water of crystallisation — caution with moisture-sensitive drug substances"),
    "Calcium_Carbonate_Precipitated": ("Calcium carbonate, precipitated",
        "Alkaline — accelerates degradation of acid-labile drug substances"),
    "Magnesium_Carbonate": ("Magnesium carbonate", "Alkaline diluent — local pH rise"),
    "Magnesium_Oxide": ("Magnesium oxide",
        "Strongly alkaline — promotes base-catalysed hydrolysis"),
    "Powdered_Cellulose": ("Cellulose, powdered",
        "Moisture sorption accelerates solid-state hydrolysis"),
    "Silicified_MCC": ("Silicified microcrystalline cellulose",
        "MCC with colloidal silica — moisture uptake and surface adsorption"),

    # --- Lubricant / Glidant ---
    "Calcium_Stearate": ("Calcium stearate",
        "Free stearic acid — amide formation with amine drug substances; reduces dissolution"),
    "Zinc_Stearate": ("Zinc stearate", "Zn²⁺ — chelation and oxidation catalysis"),
    "Sodium_Stearate": ("Sodium stearate",
        "Alkaline soap — promotes base-catalysed hydrolysis"),
    "Glyceryl_Behenate": ("Glyceryl dibehenate (Compritol)",
        "Fatty acid ester — transesterification and hydrolysis"),
    "Sodium_Benzoate_Lubricant": ("Sodium benzoate (lubricant grade)",
        "Liberates benzoic acid at acidic pH — amide formation with amines"),

    # --- Disintegrant ---
    "Alginic_Acid_Disintegrant": ("Alginic acid (disintegrant grade)",
        "Polysaccharide acid — lowers pH and exchanges ions"),
    "Polacrilin_Potassium": ("Polacrilin potassium (Amberlite IRP-88)",
        "Ion-exchange resin — ion pairing with basic drug substances; delayed dissolution"),

    # --- Polymer / Film coating ---
    "Polyethylene_Glycol_400": ("Macrogol 400 (PEG 400)",
        "🔥 Formic acid, formaldehyde and peroxides — formylation and oxidation"),
    "Polyethylene_Glycol_6000": ("Macrogol 6000 (PEG 6000)",
        "🔥 Residual formic acid and peroxides — NDSRI precursor"),
    "Povidone_K12": ("Povidone K12",
        "🔥 Residual peroxides — promotes oxidative impurities (low molecular weight grade)"),
    "Povidone_K17": ("Povidone K17", "🔥 Residual peroxides — promotes oxidative impurities"),
    "Povidone_K25": ("Povidone K25", "🔥 Residual peroxides — promotes oxidative impurities"),
    "Copovidone": ("Copovidone (Kollidon VA64)",
        "🔥 Vinylpyrrolidone/vinyl acetate copolymer — peroxides plus free acetic acid"),
    "Polyvinyl_Alcohol": ("Poly(vinyl alcohol)",
        "Residual acetate and aldehyde — reacts with amine drug substances"),
    "Methylcellulose": ("Methylcellulose",
        "Retains moisture — accelerates solid-state hydrolysis"),
    "Hydroxyethyl_Cellulose": ("Hydroxyethylcellulose",
        "Hydrophilic polymer — raises water activity"),
    "Eudragit_E": ("Aminoalkyl methacrylate copolymer E (Eudragit E)",
        "Contains a tertiary amine — N-oxide formation in the presence of peroxide"),
    "Eudragit_RS": ("Ammonio methacrylate copolymer (Eudragit RS)",
        "Quaternary ammonium — ion pairing with anionic drug substances"),
    "Shellac": ("Shellac",
        "Natural resin acid — free carboxylic acid reacts with amines; polymerises on storage"),
    "Carbomer": ("Carbomer (Carbopol)",
        "Poly(acrylic acid) — strongly acidic; ion pairing with basic drug substances"),

    # --- Surfactant ---
    "Polysorbate_40": ("Polysorbate 40 (Tween 40)",
        "Ester surfactant — peroxide formation and hydrolysis"),
    "Sorbitan_Monooleate": ("Sorbitan monooleate (Span 80)",
        "🔥 Unsaturated fatty acid ester — autoxidation peroxides"),
    "Docusate_Sodium": ("Docusate sodium",
        "Sulfosuccinate ester — hydrolysis liberates the free acid"),
    "Cetrimide": ("Cetrimide",
        "Quaternary ammonium — ion-pair impurities with anionic drug substances"),
    "Poloxamer_188": ("Poloxamer 188 (Pluronic F68)",
        "Polyether — autoxidation peroxides and formic acid"),

    # --- Buffer / pH modifier ---
    "Fumaric_Acid": ("Fumaric acid",
        "🔥 α,β-Unsaturated diacid — Michael addition with amine and thiol drug substances"),
    "Maleic_Acid": ("Maleic acid",
        "🔥 α,β-Unsaturated diacid — Michael addition; strongly acidic"),
    "Succinic_Acid": ("Succinic acid",
        "Diacid — salt formation and amide condensation with amines"),
    "Malic_Acid": ("Malic acid", "Hydroxy acid — esterification and amide formation"),
    "Lactic_Acid": ("Lactic acid", "Hydroxy acid — lactide formation and esterification"),
    "Acetic_Acid_Glacial": ("Acetic acid, glacial",
        "Acetylates amine drug substances; esterifies alcohols"),
    "Sodium_Acetate": ("Sodium acetate",
        "Buffer salt — free acetic acid can acetylate amines"),
    "Meglumine": ("Meglumine",
        "🔥 Secondary amine with a polyol chain — nitrosation risk and Maillard reaction"),
    "Arginine": ("Arginine",
        "Primary amine and guanidine — Maillard reaction and amide formation"),
    "Glycine": ("Glycine", "Primary amine — Maillard reaction and amino acid adducts"),

    # --- Sweetener / Flavour ---
    "Acesulfame_Potassium": ("Acesulfame potassium",
        "Sulfamate derivative — hydrolysis under acidic conditions"),
    "Neotame": ("Neotame",
        "Peptide ester — hydrolysis and diketopiperazine formation"),
    "Sodium_Cyclamate": ("Sodium cyclamate",
        "Sulfamate — caution with nitrite (cyclohexylamine liberation)"),
    "Vanillin": ("Vanillin",
        "🔥 Aromatic aldehyde — Schiff base with primary amine drug substances"),
    "Menthol": ("Menthol", "Terpene alcohol — eutectic mixture formation"),

    # --- Preservative / Antioxidant ---
    "Butylparaben": ("Butyl parahydroxybenzoate",
        "Ester preservative — aminolysis and hydrolysis"),
    "Ethylparaben": ("Ethyl parahydroxybenzoate",
        "Ester preservative — aminolysis and hydrolysis"),
    "Sorbic_Acid": ("Sorbic acid",
        "🔥 Conjugated unsaturated acid — Michael addition and oxidative browning"),
    "Sodium_Sulfite": ("Sodium sulfite",
        "🔥 Reducing agent — bisulfite adducts with aldehyde and ketone drug substances"),
    "Sodium_Bisulfite": ("Sodium bisulfite",
        "🔥 Reducing agent — carbonyl adducts and disulfide reduction"),
    "Sodium_Thiosulfate": ("Sodium thiosulfate",
        "Reducing agent — reduces oxidising drug substances and cleaves disulfides"),
    "Propyl_Gallate": ("Propyl gallate",
        "Polyphenol ester — metal complexation and oxidative discolouration"),
    "Alpha_Tocopherol": ("α-Tocopherol (vitamin E)",
        "Lipophilic antioxidant — quinone by-products on oxidation"),
    "Thimerosal": ("Thiomersal",
        "Organomercurial — thiol exchange with thiol drug substances; restricted use"),
    "Chlorocresol": ("Chlorocresol",
        "Halogenated phenol — polymer sorption and oxidative discolouration"),

    # --- Solvent / Co-solvent ---
    "Polyethylene_Glycol_300": ("Macrogol 300 (PEG 300)",
        "🔥 Formic acid and peroxides — formylation and oxidation"),
    "Dimethyl_Sulfoxide": ("Dimethyl sulfoxide (DMSO)",
        "May liberate formaldehyde on heating under acidic conditions"),
    "N_Methylpyrrolidone": ("N-Methylpyrrolidone (NMP)",
        "Lactam solvent — hydrolyses to methylaminobutyric acid"),
    "Benzyl_Alcohol_Solvent": ("Benzyl alcohol (solvent grade)",
        "Oxidises to benzaldehyde — Schiff base with primary amines"),
    "Castor_Oil_Polyoxyl_35": ("Macrogolglycerol ricinoleate (Cremophor EL)",
        "🔥 Unsaturated fatty acid ester — peroxide accumulation"),
}


# 카탈로그 표 컬럼 — 엔진이 한글 키로 돌려줍니다
COLUMN_EN: dict[str, str] = {
    "부형제": "Excipient",
    "카테고리": "Category",
    "출처": "Source",
    "주의보": "Caution",
    "부형제키": "Excipient",
    "Nitrite_리스크": "Nitrite risk",
    "잔류_ppm": "Residual level (ppm)",
    "가이드": "Note",
}

# 아질산염 프로파일 안내 문구
NITRITE_NOTE_EN: dict[str, str] = {
    "Crospovidone": "🔥 Povidone polymerisation and storage accumulate trace nitrite and peroxide",
    "Povidone_K30": "Residual nitrite from polymerisation; peroxide accumulation on storage",
    "Povidone_K90": "Residual nitrite from polymerisation; peroxide accumulation on storage",
    "Polyethylene_Glycol_4000": "Macrogol oxidation products may carry trace nitrite",
    "Polyethylene_Glycol_6000": "Macrogol oxidation products may carry trace nitrite",
    "Microcrystalline_Cellulose": "Trace inorganic nitrite may persist from natural pulp processing",
    "Lactose_Monohydrate": "Negligible nitrite; Maillard pathway dominates",
    "Lactose_Anhydrous": "Negligible nitrite; Maillard pathway dominates",
    "Magnesium_Stearate": "Nitrite content negligible",
}


def name_en(key: str, fallback: str = "") -> str:
    """영문 표시 이름. 표에 없으면 키를 사람이 읽을 수 있게 바꿉니다."""
    entry = EXCIPIENT_EN.get(key)
    if entry:
        return entry[0]
    if key:
        return key.replace("_", " ")
    return fallback


def alert_en(key: str, fallback: str = "") -> str:
    entry = EXCIPIENT_EN.get(key)
    return entry[1] if entry else fallback


def column_en(name: str) -> str:
    return COLUMN_EN.get(name, name)


def nitrite_note_en(key: str, fallback: str = "") -> str:
    return NITRITE_NOTE_EN.get(key, fallback)


def coverage() -> dict[str, int]:
    """번역 적용 범위 — 검증용."""
    try:
        from p72_excipient_interaction_engine import EXCIPIENT_LIBRARY

        total = len(EXCIPIENT_LIBRARY)
        covered = sum(1 for k in EXCIPIENT_LIBRARY if k in EXCIPIENT_EN)
        return {"total": total, "covered": covered, "missing": total - covered}
    except Exception:
        return {"total": 0, "covered": len(EXCIPIENT_EN), "missing": 0}
