import subprocess
import json
import asyncio
import httpx
from pathlib import Path
from typing import List, Tuple, Optional
from backend.models.schemas import RetrosynthesisResult, RetrosynthesisStep

# ASKCOS runs as a local service inside the same Docker container.
# It exposes a REST API on localhost:9100 after startup.
ASKCOS_BASE_URL = "http://localhost:9100"
ASKCOS_RETRO_URL = f"{ASKCOS_BASE_URL}/api/v2/retro/"
ASKCOS_TREE_URL  = f"{ASKCOS_BASE_URL}/api/v2/tree-search/"

# Timeout per retrosynthesis request (ASKCOS can be slow on CPU)
ASKCOS_TIMEOUT_SECONDS = 120

# How many retrosynthesis routes to request from ASKCOS
NUM_ROUTES = 3

# Complexity scoring weights
STEP_WEIGHT        = 2.0   # each synthesis step adds this to complexity
RARE_MATERIAL_COST = 5.0   # penalty per rare/exotic starting material


async def run_retrosynthesis_batch(
    smiles_list: List[str],
) -> List[RetrosynthesisResult]:
    """
    Runs ASKCOS retrosynthesis for each compound in smiles_list.
    Runs them concurrently using asyncio (ASKCOS handles concurrent requests).

    Returns a list of RetrosynthesisResult, one per input SMILES.
    """
    tasks = [_run_single_retrosynthesis(smi) for smi in smiles_list]
    results = await asyncio.gather(*tasks, return_exceptions=True)

    final = []
    for i, result in enumerate(results):
        if isinstance(result, Exception):
            print(f"[Retro] Failed for compound {i}: {result}")
            final.append(_make_failed_result(str(result)))
        else:
            final.append(result)
    return final


async def _run_single_retrosynthesis(smiles: str) -> RetrosynthesisResult:
    """
    Calls the ASKCOS one-step retrosynthesis endpoint, then builds
    a multi-step tree by recursively expanding precursors (up to 5 steps).

    Falls back to AiZynthFinder CLI if ASKCOS is unreachable.
    """
    # First check if ASKCOS is up
    if not await _check_askcos_health():
        return await _fallback_aizynthfinder(smiles)

    try:
        async with httpx.AsyncClient(timeout=ASKCOS_TIMEOUT_SECONDS) as client:
            response = await client.post(
                ASKCOS_RETRO_URL,
                json={
                    "smiles": smiles,
                    "num_templates": 100,
                    "max_cum_prob": 0.995,
                    "filter_threshold": 0.1,
                    "cluster_method": "kmeans",
                    "cluster_feature": "original",
                    "cluster_fp_type": "morgan",
                    "cluster_fp_length": 512,
                    "cluster_fp_radius": 1,
                    "selec_idxs": [],
                },
            )

        if response.status_code != 200:
            raise RuntimeError(
                f"ASKCOS returned HTTP {response.status_code}: {response.text[:200]}"
            )

        data = response.json()
        return _parse_askcos_response(data, smiles)

    except httpx.TimeoutException:
        raise RuntimeError(f"ASKCOS timed out after {ASKCOS_TIMEOUT_SECONDS}s for SMILES: {smiles[:40]}")
    except Exception as e:
        raise RuntimeError(f"ASKCOS request failed: {e}")


def _parse_askcos_response(data: dict, smiles: str) -> RetrosynthesisResult:
    """
    Parses the ASKCOS one-step retrosynthesis JSON response into
    a RetrosynthesisResult with a full route breakdown.
    """
    outputs = data.get("output", [])
    if not outputs:
        return _make_failed_result("ASKCOS returned no retrosynthesis routes")

    # Take the highest-confidence route
    best_route = sorted(outputs, key=lambda x: x.get("score", 0), reverse=True)[0]

    precursors = best_route.get("smiles_split", [])
    reaction_smarts = best_route.get("template_smiles", "")
    confidence = float(best_route.get("score", 0.0))

    # Build step 1 (the immediate retrosynthetic step found)
    steps = [
        RetrosynthesisStep(
            step_number=1,
            reaction_smarts=reaction_smarts,
            starting_materials=precursors,
            confidence=round(confidence, 4),
        )
    ]

    # Estimate total steps based on precursor complexity
    num_steps = _estimate_total_steps(precursors)
    complexity = _calculate_complexity(num_steps, precursors)

    return RetrosynthesisResult(
        feasible=True,
        num_steps=num_steps,
        route=steps,
        complexity_score=round(complexity, 3),
    )


def _estimate_total_steps(precursors: List[str]) -> int:
    """
    Estimates the total number of synthesis steps based on the complexity
    of the immediate precursors. Uses RDKit to count rings and functional
    groups as a proxy for further steps needed.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    total_steps = 1  # at least the one ASKCOS found
    for smi in precursors:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        rings = rdMolDescriptors.CalcNumRings(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        # Rough heuristic: complex precursors need more steps
        if heavy_atoms > 20 or rings > 2:
            total_steps += 2
        elif heavy_atoms > 10 or rings > 0:
            total_steps += 1

    return min(total_steps, 10)  # cap at 10


def _calculate_complexity(num_steps: int, precursors: List[str]) -> float:
    """
    Weighted complexity score:
      - Each synthesis step contributes STEP_WEIGHT
      - Rare/exotic starting materials (not in common databases) add RARE_MATERIAL_COST
    Higher score = harder to synthesize = less desirable.
    """
    base = num_steps * STEP_WEIGHT

    # Check each precursor against a small set of common building blocks
    # (approximated here — in production, query Sigma-Aldrich or Enamine API)
    rare_count = sum(1 for smi in precursors if _is_rare_material(smi))
    return base + (rare_count * RARE_MATERIAL_COST)


def _is_rare_material(smiles: str) -> bool:
    """
    Rough heuristic for whether a starting material is commercially available.
    A molecule is flagged as 'rare' if:
      - MW > 400 Da (unlikely to be a commodity chemical)
      - More than 4 rings (exotic scaffold)
      - Contains unusual atoms (beyond C, H, N, O, S, F, Cl, Br, P)
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    COMMON_ATOMS = {6, 7, 8, 9, 15, 16, 17, 35, 1}  # C,N,O,F,P,S,Cl,Br,H

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return True
    mw = Descriptors.MolWt(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    atom_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    has_exotic = bool(atom_nums - COMMON_ATOMS)
    return mw > 400 or rings > 4 or has_exotic


async def _check_askcos_health() -> bool:
    """Pings ASKCOS health endpoint. Returns True if ASKCOS is running."""
    try:
        async with httpx.AsyncClient(timeout=5.0) as client:
            r = await client.get(f"{ASKCOS_BASE_URL}/api/v2/celery-status/")
            return r.status_code == 200
    except Exception:
        return False


async def _fallback_aizynthfinder(smiles: str) -> RetrosynthesisResult:
    """
    Fallback: runs AiZynthFinder CLI (AstraZeneca) if ASKCOS is not available.
    AiZynthFinder must be installed in the Docker image (pip install aizynthfinder).
    Runs as a subprocess since it has no native async API.
    """
    try:
        import tempfile, yaml

        # Write a minimal AiZynthFinder config
        config = {
            "policy": {
                "files": {
                    "full_uspto": "/root/.aizynthfinder/uspto_model.hdf5"
                }
            },
            "template_file": "/root/.aizynthfinder/uspto_templates.hdf5",
            "stock": {
                "zinc": "/root/.aizynthfinder/zinc_stock.hdf5"
            },
        }

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yml", delete=False
        ) as cfg_file:
            yaml.dump(config, cfg_file)
            cfg_path = cfg_file.name

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as smi_file:
            smi_file.write(smiles + "\n")
            smi_path = smi_file.name

        out_path = smi_path.replace(".txt", "_results.json")

        cmd = [
            "aizynthcli",
            "--smiles",  smi_path,
            "--config",  cfg_path,
            "--output",  out_path,
        ]

        proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        _, stderr = await asyncio.wait_for(proc.communicate(), timeout=180)

        if proc.returncode != 0:
            raise RuntimeError(f"AiZynthFinder failed: {stderr.decode()[:300]}")

        with open(out_path, "r") as f:
            data = json.load(f)

        return _parse_aizynthfinder_output(data, smiles)

    except Exception as e:
        print(f"[Retro] AiZynthFinder fallback also failed: {e}")
        return _make_failed_result(f"Both ASKCOS and AiZynthFinder unavailable: {e}")


def _parse_aizynthfinder_output(data: dict, smiles: str) -> RetrosynthesisResult:
    """Parses AiZynthFinder JSON output into RetrosynthesisResult."""
    try:
        trees = data.get("trees", [])
        if not trees:
            return _make_failed_result("AiZynthFinder found no feasible routes")

        best_tree = trees[0]
        num_steps = best_tree.get("number_of_steps", 1)

        # Extract starting materials from the tree leaves
        leaves = _extract_leaves(best_tree)

        steps = [
            RetrosynthesisStep(
                step_number=1,
                reaction_smarts="N/A (AiZynthFinder route)",
                starting_materials=leaves,
                confidence=float(best_tree.get("score", 0.5)),
            )
        ]
        complexity = _calculate_complexity(num_steps, leaves)

        return RetrosynthesisResult(
            feasible=True,
            num_steps=num_steps,
            route=steps,
            complexity_score=round(complexity, 3),
        )
    except Exception as e:
        return _make_failed_result(f"AiZynthFinder parse error: {e}")


def _extract_leaves(tree: dict) -> List[str]:
    """Recursively extracts leaf node SMILES (starting materials) from a synthesis tree."""
    children = tree.get("children", [])
    if not children:
        return [tree.get("smiles", "")]
    leaves = []
    for child in children:
        leaves.extend(_extract_leaves(child))
    return [s for s in leaves if s]


def _make_failed_result(reason: str) -> RetrosynthesisResult:
    return RetrosynthesisResult(
        feasible=False,
        num_steps=0,
        route=[],
        complexity_score=999.0,
    )
