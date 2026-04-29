"""
Retrosynthesis analysis using SA Score + functional-group-driven disconnection.

Provides:
  - Synthetic feasibility (SA score, num_steps estimate)
  - BRICS fragment decomposition for starting material identification
  - Functional group detection and priority-ranked disconnections
  - Named reactions with reagents, conditions, and step-by-step protocol
  - Reaction schematics (A + B → Product)
  - Synthesis narrative overview
  
All analysis is purely RDKit-based — no external API required.
"""

import asyncio
from typing import List, Dict, Optional, Tuple, Any
from backend.models.schemas import RetrosynthesisResult, RetrosynthesisStep


# ── Reaction knowledge base ────────────────────────────────────────────────────
#
# Each entry maps a detectable functional group (SMARTS) to the
# retrosynthetic disconnection that most likely formed it.
# Priority = order in list (earlier = checked first = higher retrosynthetic value).

_REACTION_KB: List[Dict[str, Any]] = [
    {
        "smarts": "[NX3;H1,H2][CX3](=[OX1])",
        "fg_name": "Amide",
        "reaction_name": "Amide Coupling",
        "reaction_type": "C–N Bond Formation",
        "schematic_template": "R-COOH  +  H₂N-R'  →  R-C(=O)-NH-R'",
        "reagents": ["HATU (1.1 eq)", "DIPEA (2.0 eq)", "DMF, anhydrous", "or EDC·HCl (1.2 eq) + HOBt (1.1 eq)"],
        "conditions": "RT, dry DMF, 12–24 h, under N₂",
        "protocol_steps": [
            "Dissolve carboxylic acid (1.0 eq) in anhydrous DMF (0.1 M) under N₂.",
            "Add HATU (1.1 eq) and DIPEA (2.0 eq). Stir 10 min at RT to pre-activate.",
            "Add amine component (1.0 eq) in DMF dropwise. Stir 12–24 h at RT.",
            "Dilute with EtOAc, wash with sat. NaHCO₃ (×2), water (×1), brine (×1).",
            "Dry over MgSO₄, filter, concentrate in vacuo.",
            "Purify by flash column chromatography (SiO₂, EtOAc/hexane gradient).",
            "Confirm by ¹H NMR: amide NH appears ~8–9 ppm (broad).",
        ],
        "estimated_yield": "60–90%",
        "notes": (
            "HATU preferred for hindered or electron-poor amines. "
            "For peptide-like substrates, add HOBt to suppress epimerization. "
            "Thionyl chloride activation → acyl chloride is faster but harsher."
        ),
    },
    {
        "smarts": "[NX3;H1,H0][SX4](=[OX1])(=[OX1])",
        "fg_name": "Sulfonamide",
        "reaction_name": "Sulfonamide Formation",
        "reaction_type": "N–S Bond Formation",
        "schematic_template": "R-SO₂Cl  +  H₂N-R'  →  R-SO₂-NH-R'",
        "reagents": ["Sulfonyl chloride (1.1 eq)", "Et₃N (2.0 eq) or pyridine (5 eq)", "DCM, dry"],
        "conditions": "0°C → RT, DCM, 2–4 h",
        "protocol_steps": [
            "Cool amine (1.0 eq) in dry DCM (0.2 M) to 0°C under N₂.",
            "Add Et₃N (2.0 eq) dropwise, then sulfonyl chloride (1.1 eq) dropwise.",
            "Allow to warm to RT slowly, stir 2–4 h.",
            "Wash with 1 M HCl (×2), sat. NaHCO₃ (×1), brine (×1).",
            "Dry (MgSO₄), filter, concentrate.",
            "Purify by column (SiO₂) or recrystallize from EtOH/H₂O.",
        ],
        "estimated_yield": "70–95%",
        "notes": (
            "Pyridine as solvent/base works excellently for aromatic amines. "
            "For primary amines, N-monosubstitution can be challenging — "
            "consider protecting one N-H first."
        ),
    },
    {
        "smarts": "[OX2][CX3](=[OX1])[#6]",
        "fg_name": "Ester",
        "reaction_name": "DCC/DMAP Esterification",
        "reaction_type": "C–O Bond Formation",
        "schematic_template": "R-COOH  +  HO-R'  →  R-C(=O)-O-R'",
        "reagents": ["DCC (1.1 eq)", "DMAP (0.1 eq, cat.)", "DCM", "or H₂SO₄ (cat.) in toluene for Fischer"],
        "conditions": "0°C → RT, DCM, 12 h (DCC); or reflux/Dean-Stark (Fischer esterification)",
        "protocol_steps": [
            "Dissolve carboxylic acid (1.0 eq) and alcohol (1.05 eq) in DCM (0.1 M), cool to 0°C.",
            "Add DMAP (0.1 eq) and DCC (1.1 eq).",
            "Warm to RT, stir 12 h under N₂.",
            "Filter off precipitated DCU (dicyclohexylurea) through Celite.",
            "Wash filtrate with sat. NaHCO₃ (×2), brine.",
            "Dry (Na₂SO₄), concentrate, purify by column (SiO₂, EtOAc/hexane).",
            "Confirm ester formation: IR ~1735 cm⁻¹ C=O stretch.",
        ],
        "estimated_yield": "65–90%",
        "notes": (
            "Acyl chloride + alcohol in pyridine gives faster, higher-yielding ester. "
            "Fischer esterification (acid + alcohol, H₂SO₄, reflux) works for simple esters."
        ),
    },
    {
        "smarts": "[OX2H][CX4;H1]([#6])[#6]",
        "fg_name": "Secondary Alcohol",
        "reaction_name": "Ketone Reduction (NaBH₄)",
        "reaction_type": "Reduction (C=O → CHOH)",
        "schematic_template": "R-C(=O)-R'  +  NaBH₄  →  R-CH(OH)-R'",
        "reagents": ["NaBH₄ (2.0 eq)", "MeOH or EtOH"],
        "conditions": "0°C → RT, MeOH, 1–2 h",
        "protocol_steps": [
            "Dissolve ketone (1.0 eq) in MeOH (0.2 M), cool to 0°C.",
            "Add NaBH₄ (2.0 eq) portionwise over 10 min.",
            "Stir at 0°C for 30 min, then warm to RT for 1 h.",
            "Quench carefully with sat. NH₄Cl (dropwise — exothermic).",
            "Extract with EtOAc (×3), wash brine, dry (MgSO₄), concentrate.",
            "Purify by column or use crude if pure enough by TLC.",
        ],
        "estimated_yield": "80–98%",
        "notes": (
            "Highly reliable. LiAlH₄ in THF for less reactive ketones or ester reduction. "
            "CBS catalyst for asymmetric reduction if stereochemistry is critical."
        ),
    },
    {
        "smarts": "[OX2;!$(OC=O);!$(Oc)]([CX4])[CX4]",
        "fg_name": "Aliphatic Ether",
        "reaction_name": "Williamson Ether Synthesis",
        "reaction_type": "C–O Bond Formation (SN2)",
        "schematic_template": "R-OH  +  R'-X  →  R-O-R'",
        "reagents": ["NaH (1.2 eq) or K₂CO₃ (2.0 eq)", "Alkyl halide (1.1 eq)", "DMF (NaH) or acetone (K₂CO₃)"],
        "conditions": "0°C → RT for NaH/DMF; 60°C for K₂CO₃/acetone, 4–12 h",
        "protocol_steps": [
            "Suspend NaH (1.2 eq, 60% in oil — wash with dry hexane ×2 first) in dry DMF under N₂, cool to 0°C.",
            "Add alcohol (1.0 eq) in DMF dropwise. Stir 30 min until H₂ evolution ceases.",
            "Add alkyl halide (1.1 eq) dropwise at 0°C.",
            "Warm to RT, stir 4–12 h (monitor by TLC).",
            "Quench carefully with sat. NH₄Cl. Extract with EtOAc (×3).",
            "Wash with water (×2), brine. Dry (MgSO₄), concentrate, purify.",
        ],
        "estimated_yield": "50–85%",
        "notes": (
            "NaH is pyrophoric — handle with care. "
            "K₂CO₃ in acetone/DMF is safer and works well for benzyl and allyl halides. "
            "Mitsunobu reaction (PPh₃, DIAD) for inversion of alcohol stereochemistry."
        ),
    },
    {
        "smarts": "[NX3;H0;!$(NC=O);!$(NS=O)]([CX4])([CX4])[CX4]",
        "fg_name": "Tertiary Amine",
        "reaction_name": "Reductive Amination",
        "reaction_type": "C–N Bond Formation",
        "schematic_template": "R₂NH  +  R'-CHO  →  R₂N-CH₂-R'  (via iminium, then reduction)",
        "reagents": ["Aldehyde (1.1 eq)", "NaBH₃CN (1.5 eq) or NaBH(OAc)₃ (1.5 eq)", "MeOH or DCE", "AcOH (cat., 0.1 eq)"],
        "conditions": "RT, MeOH or DCE, 12–24 h",
        "protocol_steps": [
            "Combine amine (1.0 eq) and carbonyl compound (1.1 eq) in MeOH (0.2 M).",
            "Add AcOH (0.1 eq). Stir 1 h at RT to form iminium ion.",
            "Add NaBH₃CN (1.5 eq) portionwise. Stir 12 h at RT.",
            "Basify to pH ~10 with 2 M NaOH. Extract with EtOAc (×3).",
            "Wash brine, dry (Na₂SO₄), concentrate.",
            "Purify by column (SiO₂, amine modifier 1% Et₃N recommended).",
        ],
        "estimated_yield": "55–85%",
        "notes": (
            "NaBH(OAc)₃ in DCE is milder and more chemoselective. "
            "Avoid ketones with bulky groups — iminium formation is slow. "
            "Leuckart reaction for formylation with formic acid is an alternative."
        ),
    },
    {
        "smarts": "[NX3][CX3](=[OX1])[NX3]",
        "fg_name": "Urea",
        "reaction_name": "Urea Formation (Isocyanate Route)",
        "reaction_type": "C–N Bond Formation",
        "schematic_template": "R-NH₂  +  R'-NCO  →  R-NH-C(=O)-NH-R'",
        "reagents": ["Isocyanate (1.0 eq)", "Et₃N (1.1 eq, optional)", "DCM or THF"],
        "conditions": "RT, DCM, 1–4 h",
        "protocol_steps": [
            "Dissolve amine (1.0 eq) in DCM (0.2 M) under N₂.",
            "Add isocyanate (1.0 eq) dropwise at RT (exothermic for reactive amines).",
            "Stir 1–4 h. Urea often precipitates from solution.",
            "If solid: filter, wash with cold DCM, dry under vacuum.",
            "If in solution: wash with water, dry (MgSO₄), concentrate, recrystallize.",
        ],
        "estimated_yield": "70–95%",
        "notes": (
            "CDI (carbonyldiimidazole) + two different amines gives unsymmetrical ureas cleanly. "
            "Phosgene route also works but requires extreme caution (toxic gas). "
            "Ureas are often low-solubility solids — good for crystallization-based purification."
        ),
    },
    {
        "smarts": "[NX3][CX3](=[OX1])[OX2]",
        "fg_name": "Carbamate",
        "reaction_name": "Boc Protection / Carbamate Synthesis",
        "reaction_type": "N–CO Bond Formation",
        "schematic_template": "R-NH₂  +  Boc₂O  →  R-NH-C(=O)-O-tBu",
        "reagents": ["Boc₂O (1.1 eq)", "Et₃N (1.5 eq) or DMAP (cat.)", "DCM"],
        "conditions": "0°C → RT, DCM, 2–4 h",
        "protocol_steps": [
            "Dissolve amine (1.0 eq) in DCM (0.2 M), cool to 0°C.",
            "Add Et₃N (1.5 eq) dropwise, then Boc₂O (1.1 eq) in DCM dropwise.",
            "Warm to RT, stir 2–4 h (monitor by TLC).",
            "Wash with 10% citric acid (×2), sat. NaHCO₃ (×1), brine.",
            "Dry (MgSO₄), filter, concentrate.",
            "Boc removal: TFA/DCM 1:1, 30 min RT, then evaporate TFA completely.",
        ],
        "estimated_yield": "75–95%",
        "notes": (
            "Boc is acid-labile (TFA/DCM). Cbz is H₂/Pd-labile. "
            "Fmoc is base-labile (piperidine 20% in DMF). "
            "Choose protecting group based on downstream transformations."
        ),
    },
    {
        "smarts": "[NX3;H1,H2][c]",
        "fg_name": "Aryl Amine (Aniline)",
        "reaction_name": "Buchwald-Hartwig Amination",
        "reaction_type": "Pd-Catalysed C(Ar)–N Cross-Coupling",
        "schematic_template": "Ar-X  +  HN-R₂  →  Ar-N-R₂  (Pd₂dba₃/XPhos catalyst)",
        "reagents": ["Pd₂(dba)₃ (2 mol%)", "XPhos or BINAP (4 mol%)", "Cs₂CO₃ (2.0 eq)", "Toluene or 1,4-dioxane, dry"],
        "conditions": "100–110°C, toluene, 12–24 h, under N₂",
        "protocol_steps": [
            "Charge aryl halide (1.0 eq), amine (1.2 eq), Cs₂CO₃ (2.0 eq) in a sealed vial with toluene (0.2 M).",
            "Degas thoroughly: three freeze-pump-thaw cycles or N₂ sparge 20 min.",
            "Add Pd₂(dba)₃ (2 mol%) and XPhos ligand (4 mol%).",
            "Seal vial, heat to 100–110°C with stirring, 12–24 h.",
            "Cool to RT, dilute with EtOAc, filter through Celite.",
            "Wash organic layer with water (×2), brine. Dry (MgSO₄), concentrate.",
            "Purify by column (SiO₂, EtOAc/hexane). Confirm by LC-MS.",
        ],
        "estimated_yield": "50–85%",
        "notes": (
            "SNAr works without Pd if the aryl ring carries strong EWG (NO₂, CN) ortho/para to halide. "
            "Aryl bromides and iodides are preferred substrates; chlorides need bulky ligands (SPhos). "
            "RuPhos-Pd-G3 precatalyst is an excellent, bench-stable alternative."
        ),
    },
    {
        "smarts": "[cX3](-[cX3]1ccccc1)-[cX3]1ccccc1",
        "fg_name": "Biaryl",
        "reaction_name": "Suzuki-Miyaura Coupling",
        "reaction_type": "Pd-Catalysed C(Ar)–C(Ar) Cross-Coupling",
        "schematic_template": "Ar-X  +  Ar'-B(OH)₂  →  Ar-Ar'  (Pd(PPh₃)₄ catalyst)",
        "reagents": ["Pd(PPh₃)₄ (3 mol%) or Pd(dppf)Cl₂·DCM (2 mol%)", "K₂CO₃ or Na₂CO₃ (2.0 eq)", "DMF/H₂O 4:1 or toluene/EtOH"],
        "conditions": "80–90°C, DMF/H₂O (4:1), 4–12 h, under N₂",
        "protocol_steps": [
            "Combine aryl halide (1.0 eq) and arylboronic acid (1.2 eq) in DMF/H₂O 4:1 (0.2 M).",
            "Add K₂CO₃ (2.0 eq). Sparge with N₂ for 15 min.",
            "Add Pd(PPh₃)₄ (3 mol%) under N₂ counter-flow.",
            "Heat to 80–90°C, stir 4–12 h under N₂ balloon.",
            "Cool, dilute with EtOAc. Wash water (×3), brine. Dry (MgSO₄), concentrate.",
            "Purify by column chromatography (SiO₂).",
        ],
        "estimated_yield": "60–90%",
        "notes": (
            "Aryl bromides and iodides react readily at 80°C. "
            "Aryl chlorides need bulky phosphine ligands (SPhos, RuPhos) or at higher temps. "
            "MIDA boronates are air-stable alternatives to boronic acids."
        ),
    },
    {
        "smarts": "[CX2]#[NX1]",
        "fg_name": "Nitrile",
        "reaction_name": "Pd-Catalysed Cyanation",
        "reaction_type": "C–CN Bond Formation",
        "schematic_template": "Ar-X  +  Zn(CN)₂  →  Ar-CN  (Pd₂dba₃/dppf catalyst)",
        "reagents": ["Zn(CN)₂ (0.6 eq)", "Pd₂(dba)₃ (2 mol%)", "dppf or dppb (4 mol%)", "DMF, dry"],
        "conditions": "100°C, DMF, 12 h, N₂",
        "protocol_steps": [
            "Combine aryl halide (1.0 eq) and Zn(CN)₂ (0.6 eq) in dry DMF (0.2 M).",
            "Add Pd₂(dba)₃ (2 mol%) and dppf ligand (4 mol%). Sparge N₂ 15 min.",
            "Heat to 100°C, stir 12 h under N₂.",
            "Cool, dilute EtOAc, wash sat. NH₄Cl (×3), brine. Dry, concentrate.",
            "Purify by column (SiO₂, EtOAc/hexane).",
        ],
        "estimated_yield": "40–80%",
        "notes": (
            "Rosenmund-von Braun (CuCN, NMP, 170°C) works without Pd but requires harsh conditions. "
            "Pd route preferred for sensitive substrates. "
            "Nitrile can be hydrolysed → amide → amine or carboxylic acid downstream."
        ),
    },
    {
        "smarts": "[#6][FX1,ClX1,BrX1,IX1]",
        "fg_name": "Alkyl/Aryl Halide",
        "reaction_name": "Halogenation (Radical or Electrophilic)",
        "reaction_type": "C–X Bond Formation",
        "schematic_template": "R-H  +  X₂  →  R-X  (radical or EAS)",
        "reagents": ["NBS (1.05 eq) for benzylic/allylic", "Br₂ (1.0 eq) for aromatic EAS", "hν or AIBN initiator"],
        "conditions": "CCl₄ or DCM, hν or 80°C (radical); 0°C, DCM (EAS)",
        "protocol_steps": [
            "Dissolve substrate (1.0 eq) in CCl₄ (0.2 M).",
            "Add NBS (1.05 eq) and a catalytic amount of AIBN or benzoyl peroxide.",
            "Irradiate with a 500W lamp or heat to 80°C, stir 2–4 h (monitor by TLC).",
            "Cool, filter off succinimide byproduct.",
            "Wash filtrate with sat. NaHCO₃, brine. Dry (MgSO₄), concentrate.",
            "Purify by column or recrystallize.",
        ],
        "estimated_yield": "50–80%",
        "notes": (
            "NBS is the standard for benzylic bromination. "
            "Aromatic halogenation (EAS): use Br₂ + FeBr₃ or AlCl₃ as Lewis acid. "
            "Halide is a key coupling partner for Pd cross-couplings downstream."
        ),
    },
]


# ── Async entry point ─────────────────────────────────────────────────────────

async def run_retrosynthesis_batch(
    smiles_list: List[str],
) -> List[RetrosynthesisResult]:
    """
    SA Score-based feasibility + functional-group disconnection analysis.
    Pure RDKit — no external service required.
    """
    results = []
    for i, smi in enumerate(smiles_list):
        try:
            result = _assess_single(smi)
            results.append(result)
        except Exception as e:
            print(f"[Retro] Failed for compound {i}: {e}")
            results.append(_make_failed_result(str(e)))
    return results


# ── Core analysis ─────────────────────────────────────────────────────────────

def _assess_single(smiles: str) -> RetrosynthesisResult:
    from rdkit import Chem, RDConfig
    from rdkit.Chem import rdMolDescriptors
    import os, sys

    # Load SA Score contrib module (ships with all RDKit distributions)
    sa_score_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if sa_score_path not in sys.path:
        sys.path.append(sa_score_path)
    import sascorer

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return _make_failed_result(f"Invalid SMILES: {smiles}")

    sa_score = sascorer.calculateScore(mol)
    print(f"[Retro] SA Score for {smiles[:50]}: {sa_score:.3f}")

    # ── Feasibility thresholds (literature standard) ──────────────────────
    # 1.0 – 3.5 : Easy (commercial building blocks, 1-3 steps)
    # 3.5 – 5.0 : Moderate (4-6 steps, needs skilled chemist)
    # 5.0 – 6.0 : Hard (6-8 steps, multi-step campaign)
    # 6.0 – 10.0: Very hard / practically infeasible

    feasible = sa_score <= 6.5

    if sa_score <= 2.0:
        num_steps, category = 1, "easy"
    elif sa_score <= 3.5:
        num_steps, category = 3, "easy"
    elif sa_score <= 5.0:
        num_steps, category = 5, "moderate"
    elif sa_score <= 6.0:
        num_steps, category = 7, "hard"
    elif sa_score <= 6.5:
        num_steps, category = 9, "very hard"
    else:
        num_steps, category = 0, "infeasible"

    # ── Detect functional groups ──────────────────────────────────────────
    detected_fgs = _detect_functional_groups(mol)

    # ── Build disconnection-based synthesis steps ─────────────────────────
    starting_materials = _fragment_molecule(mol)
    canonical_smiles = Chem.MolToSmiles(mol)

    steps = _build_synthesis_steps(
        mol=mol,
        detected_fgs=detected_fgs,
        product_smiles=canonical_smiles,
        starting_materials=starting_materials,
        sa_score=sa_score,
        feasible=feasible,
    )

    # ── Build human-readable synthesis overview ───────────────────────────
    overview = _build_overview(
        smiles=canonical_smiles,
        sa_score=sa_score,
        category=category,
        num_steps=num_steps,
        detected_fgs=detected_fgs,
        starting_materials=starting_materials,
    )

    # ── Alternative routes hint ───────────────────────────────────────────
    alt_hint = _alt_routes_hint(detected_fgs)

    # Complexity score 0–100 for pipeline scorer
    complexity = round(sa_score * 10.0, 1)

    return RetrosynthesisResult(
        feasible=feasible,
        num_steps=num_steps,
        route=steps,
        complexity_score=complexity,
        sa_score=round(sa_score, 3),
        complexity_category=category,
        functional_groups_detected=[fg["fg_name"] for fg in detected_fgs],
        disconnection_points=[fg["reaction_type"] for fg in detected_fgs],
        synthesis_overview=overview,
        alternative_routes_hint=alt_hint,
    )


# ── Functional group detection ────────────────────────────────────────────────

def _detect_functional_groups(mol) -> List[Dict[str, Any]]:
    """
    Match SMARTS patterns against the molecule in priority order.
    Returns up to 4 matched entries (avoids over-reporting overlapping groups).
    """
    from rdkit import Chem

    matched: List[Dict[str, Any]] = []
    seen_fg_names = set()

    for entry in _REACTION_KB:
        if entry["fg_name"] in seen_fg_names:
            continue
        try:
            pattern = Chem.MolFromSmarts(entry["smarts"])
            if pattern is None:
                continue
            if mol.HasSubstructMatch(pattern):
                matched.append(entry)
                seen_fg_names.add(entry["fg_name"])
                if len(matched) >= 4:
                    break
        except Exception:
            continue

    return matched


# ── Synthesis step builder ────────────────────────────────────────────────────

def _build_synthesis_steps(
    mol,
    detected_fgs: List[Dict[str, Any]],
    product_smiles: str,
    starting_materials: List[str],
    sa_score: float,
    feasible: bool,
) -> List[RetrosynthesisStep]:
    """
    Build a list of RetrosynthesisStep objects from detected functional groups.
    Each step covers one key bond disconnection.
    """
    from rdkit import Chem

    if not feasible:
        return []

    steps: List[RetrosynthesisStep] = []
    confidence_base = round(1.0 - (sa_score / 10.0), 3)
    total_fgs = len(detected_fgs)

    for idx, fg in enumerate(detected_fgs):
        step_num = idx + 1

        # Schematic: use template with generic labels
        schematic = fg.get("schematic_template", "A + B → Product")

        # Starting materials: use BRICS fragments for step 1, generic labels after
        if idx == 0:
            sm_list = starting_materials[:3] if starting_materials else ["Fragment A", "Fragment B"]
        else:
            sm_list = [f"Intermediate from Step {idx}", "Reagent"]

        # Confidence degrades slightly for later steps (more uncertainty)
        step_confidence = round(max(0.1, confidence_base - idx * 0.05), 3)

        step = RetrosynthesisStep(
            step_number=step_num,
            reaction_name=fg["reaction_name"],
            reaction_type=fg["reaction_type"],
            reaction_smarts=fg["smarts"],
            starting_materials=sm_list,
            product_smiles=product_smiles if idx == total_fgs - 1 else None,
            schematic=schematic,
            reagents=fg.get("reagents", []),
            conditions=fg.get("conditions", ""),
            protocol_steps=fg.get("protocol_steps", []),
            confidence=step_confidence,
            estimated_yield=fg.get("estimated_yield"),
            notes=fg.get("notes"),
        )
        steps.append(step)

    # Always add a final assembly/purification step
    if steps:
        steps.append(
            RetrosynthesisStep(
                step_number=len(steps) + 1,
                reaction_name="Final Purification & Characterisation",
                reaction_type="Work-up",
                reaction_smarts="",
                starting_materials=["Crude product from final coupling step"],
                product_smiles=product_smiles,
                schematic=f"Crude → Purified target: {product_smiles[:60]}",
                reagents=["SiO₂ gel (flash column)", "HPLC (reverse-phase, C18)", "rotary evaporator"],
                conditions="Flash column: gradient EtOAc/hexane. HPLC: MeCN/H₂O + 0.1% TFA.",
                protocol_steps=[
                    "Perform flash column chromatography to remove major impurities.",
                    "Check purity by HPLC (>95% required for biological testing).",
                    "Confirm identity by ¹H NMR, ¹³C NMR, and high-resolution MS (HRMS).",
                    "Determine exact mass by HRMS (ESI+): [M+H]⁺ within 5 ppm.",
                    "Record melting point if solid. Measure optical rotation if chiral.",
                    "Determine aqueous solubility (nephelometry or DMSO stock dilution assay).",
                ],
                confidence=round(confidence_base, 3),
                estimated_yield="Purification recovery: 70–95%",
                notes=(
                    "Always run a purity check before bioassay submission. "
                    "DMSO stock (10 mM) stored at -20°C is standard for assay delivery. "
                    "Document all spectral data and retain a reference sample."
                ),
            )
        )

    return steps


# ── Narrative overview ────────────────────────────────────────────────────────

def _build_overview(
    smiles: str,
    sa_score: float,
    category: str,
    num_steps: int,
    detected_fgs: List[Dict[str, Any]],
    starting_materials: List[str],
) -> str:
    fg_names = [fg["fg_name"] for fg in detected_fgs]
    rxn_types = [fg["reaction_name"] for fg in detected_fgs]

    fg_str = ", ".join(fg_names) if fg_names else "no complex functional groups"
    rxn_str = "; ".join(rxn_types) if rxn_types else "straightforward chemistry"

    sm_str = (
        ", ".join(starting_materials[:3]) + (" ..." if len(starting_materials) > 3 else "")
        if starting_materials else "commercial building blocks"
    )

    step_word = "step" if num_steps == 1 else "steps"

    feasibility_note = {
        "easy": "This compound should be readily accessible in a typical medicinal chemistry lab.",
        "moderate": "Synthesis is achievable but requires careful planning and execution.",
        "hard": "Multi-step synthesis requiring significant resource investment. Consider analogue simplification.",
        "very hard": "Highly complex; consider whether simpler analogues retain activity.",
        "infeasible": "SA Score indicates this molecule is practically infeasible to synthesise with current methods.",
    }.get(category, "")

    return (
        f"SA Score: {sa_score:.2f} / 10 — complexity category: {category.upper()} "
        f"(estimated {num_steps} synthetic {step_word}). "
        f"Key functional groups detected: {fg_str}. "
        f"Likely key transformations: {rxn_str}. "
        f"Suggested starting materials (BRICS fragments): {sm_str}. "
        f"{feasibility_note}"
    )


def _alt_routes_hint(detected_fgs: List[Dict[str, Any]]) -> Optional[str]:
    hints = {
        "Amide": "Alternative amide routes: acyl chloride + amine (faster); Weinreb amide for ketone synthesis.",
        "Biaryl": "Alternative biaryl couplings: Negishi (ArZnX), Kumada (ArMgX), Stille (ArSnBu₃).",
        "Aryl Amine (Aniline)": "Alternative aryl C–N: SNAr (no Pd, needs EWG on ring); Chan-Lam (ArB(OH)₂ + amine, Cu cat).",
        "Ester": "Alternative ester routes: acyl chloride + alcohol; enzymatic esterification for selectivity.",
        "Nitrile": "Alternative nitrile routes: Rosenmund-von Braun (CuCN, NMP, 170°C) for aryl nitriles.",
    }
    fg_names = [fg["fg_name"] for fg in detected_fgs]
    relevant = [hints[n] for n in fg_names if n in hints]
    return " | ".join(relevant) if relevant else None


# ── BRICS fragmentation ───────────────────────────────────────────────────────

def _fragment_molecule(mol) -> List[str]:
    """
    BRICS decomposition — identifies likely synthetically accessible subunits.
    """
    from rdkit.Chem import BRICS, AllChem
    from rdkit import Chem

    try:
        frags = BRICS.BRICSDecompose(mol)
        import re
        clean = []
        for frag in frags:
            # Remove all BRICS dummy atom labels like [3*], [14*]
            cleaned = re.sub(r"\[\d+\*\]", "[H]", frag)
            m = Chem.MolFromSmiles(cleaned)
            if m and m.GetNumHeavyAtoms() >= 3:
                clean.append(Chem.MolToSmiles(m))
        return clean[:5] if clean else [Chem.MolToSmiles(mol)]
    except Exception:
        return [Chem.MolToSmiles(mol)]


# ── Error result ──────────────────────────────────────────────────────────────

def _make_failed_result(reason: str) -> RetrosynthesisResult:
    print(f"[Retro] Failed: {reason}")
    return RetrosynthesisResult(
        feasible=False,
        num_steps=0,
        route=[],
        complexity_score=999.0,
        sa_score=10.0,
        complexity_category="infeasible",
        functional_groups_detected=[],
        disconnection_points=[],
        synthesis_overview=f"Retrosynthesis failed: {reason}",
        alternative_routes_hint=None,
    )