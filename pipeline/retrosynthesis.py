# backend/pipeline/retrosynthesis.py
import asyncio
from typing import List
from backend.models.schemas import RetrosynthesisResult, RetrosynthesisStep


async def run_retrosynthesis_batch(
    smiles_list: List[str],
) -> List[RetrosynthesisResult]:
    """
    Runs SA Score-based retrosynthesis assessment for each compound.
    Uses RDKit's Synthetic Accessibility (SA) Score — no external service needed.
    SA Score range: 1 (easy) to 10 (very hard). Below 3.5 = synthetically feasible.
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


def _assess_single(smiles: str) -> RetrosynthesisResult:
    from rdkit import Chem
    from rdkit.Chem import RDConfig, Descriptors, rdMolDescriptors
    import os
    import sys

    # Load SA Score contrib module (ships with RDKit)
    sa_score_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if sa_score_path not in sys.path:
        sys.path.append(sa_score_path)

    import sascorer

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return _make_failed_result(f"Invalid SMILES: {smiles}")

    sa_score = sascorer.calculateScore(mol)
    print(f"[Retro] SA Score for {smiles[:40]}: {sa_score:.3f}")

    # SA Score thresholds (standard in literature)
    # 1.0 - 3.5: easy to synthesize
    # 3.5 - 6.0: moderate difficulty
    # 6.0 - 10.0: very hard / infeasible
    feasible = sa_score <= 6.0

    # Estimate number of steps from SA score
    if sa_score <= 2.0:
        num_steps = 1
    elif sa_score <= 3.5:
        num_steps = 2
    elif sa_score <= 5.0:
        num_steps = 4
    elif sa_score <= 6.0:
        num_steps = 6
    else:
        num_steps = 0  # infeasible

    # Identify starting materials via retrosynthetic fragmentation
    starting_materials = _fragment_molecule(mol)

    steps = []
    if feasible and num_steps > 0:
        steps = [
            RetrosynthesisStep(
                step_number=1,
                reaction_smarts="SA Score based assessment",
                starting_materials=starting_materials,
                confidence=round(1.0 - (sa_score / 10.0), 3),
            )
        ]

    # Complexity score = SA score * 10 (maps to 0-100 scale for scoring)
    complexity = round(sa_score * 10, 1)

    return RetrosynthesisResult(
        feasible=feasible,
        num_steps=num_steps,
        route=steps,
        complexity_score=complexity,
    )


def _fragment_molecule(mol) -> List[str]:
    """
    Breaks molecule into likely starting material fragments using
    RDKit's BRICS decomposition — identifies synthetically accessible subunits.
    """
    from rdkit.Chem import BRICS, AllChem
    from rdkit import Chem

    try:
        frags = BRICS.BRICSDecompose(mol)
        # Clean up BRICS dummy atoms ([3*], [4*] etc.)
        clean = []
        for frag in frags:
            cleaned = frag.replace("[3*]", "[H]").replace("[4*]", "[H]") \
                         .replace("[1*]", "[H]").replace("[2*]", "[H]") \
                         .replace("[5*]", "[H]").replace("[6*]", "[H]") \
                         .replace("[7*]", "[H]").replace("[8*]", "[H]") \
                         .replace("[9*]", "[H]").replace("[10*]", "[H]")
            m = Chem.MolFromSmiles(cleaned)
            if m and m.GetNumHeavyAtoms() >= 3:
                clean.append(Chem.MolToSmiles(m))
        return clean[:5] if clean else [Chem.MolToSmiles(mol)]
    except Exception:
        return [Chem.MolToSmiles(mol)]


def _make_failed_result(reason: str) -> RetrosynthesisResult:
    print(f"[Retro] Failed: {reason}")
    return RetrosynthesisResult(
        feasible=False,
        num_steps=0,
        route=[],
        complexity_score=999.0,
    )
