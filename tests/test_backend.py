"""
Run from cadd-pipeline/ root:
    cd cadd-pipeline
    python tests/test_backend.py
"""

import sys
import json
import asyncio
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

PASS = "✅"
FAIL = "❌"
WARN = "⚠️ "

results = []

def check(label: str, passed: bool, detail: str = ""):
    symbol = PASS if passed else FAIL
    print(f"  {symbol} {label}" + (f" — {detail}" if detail else ""))
    results.append(passed)


print("\n" + "="*55)
print("  CADD PIPELINE — BACKEND TEST SUITE")
print("="*55)


# ── 1. Schema validation ──────────────────────────────────────────────────────
print("\n[1] Pydantic Schemas")
try:
    from backend.models.schemas import (
        JobRequest, PipelineSteps, BindingSiteMode, DockingSpeed,
        LipinskiResult, ADMETResult, CompoundResult
    )

    # Valid request with PDB ID
    req = JobRequest(
        smiles="OC1=C(O)C(O)=CC(=C1)C1=C(O)C(=O)C2=C(O)C(O)=CC(O)=C2O1",
        pdb_id="1ABC",
    )
    check("JobRequest with PDB ID", req.pdb_id == "1ABC")
    check("SMILES stripped and stored", len(req.smiles) > 0)
    check("Default num_analogues=25", req.num_analogues == 25)
    check("Default docking speed=balanced", req.docking_speed == DockingSpeed.balanced)

    # Valid request with pdb_content
    req2 = JobRequest(
        smiles="C1=CC=CC=C1",
        pdb_content="ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00  0.00           N",
    )
    check("JobRequest with pdb_content", req2.pdb_content is not None)

    # Invalid: missing both pdb_id and pdb_content
    try:
        bad = JobRequest(smiles="CCO")
        check("Rejects missing PDB", False, "Should have raised ValueError")
    except Exception:
        check("Rejects missing PDB", True)

    # Invalid PDB ID format
    try:
        bad2 = JobRequest(smiles="CCO", pdb_id="TOOLONG")
        check("Rejects invalid PDB ID", False)
    except Exception:
        check("Rejects invalid PDB ID", True)

except Exception as e:
    check("Schema import", False, str(e))


# ── 2. SMILES Validator ───────────────────────────────────────────────────────
print("\n[2] SMILES Validator (RDKit)")
try:
    from backend.utils.smiles_validator import (
        validate_and_canonicalize, get_basic_properties
    )

    # Myricetin SMILES
    myricetin = "OC1=C(O)C(=CC(=C1)C1=C(O)C(=O)C2=CC(O)=CC(O)=C12)O"
    valid, canon, err = validate_and_canonicalize(myricetin)
    check("Myricetin SMILES valid", valid, canon[:30] if canon else err)

    # Invalid SMILES
    valid2, _, err2 = validate_and_canonicalize("NOT_A_SMILES_XYZ")
    check("Rejects invalid SMILES", not valid2)

    # Too small
    valid3, _, _ = validate_and_canonicalize("CC")
    check("Rejects too-small molecule", not valid3)

    # Basic properties
    props = get_basic_properties(myricetin)
    check("MW computed", props["molecular_weight"] > 100, f"{props['molecular_weight']} Da")
    check("Formula computed", "C" in props["formula"], props["formula"])

except Exception as e:
    check("SMILES validator import", False, str(e))


# ── 3. Analogue Generation ────────────────────────────────────────────────────
print("\n[3] Analogue Generation (RDKit BRICS)")
try:
    from backend.pipeline.analogue_gen import generate_analogues

    # Use a simpler molecule for faster testing
    test_smiles = "c1ccc(O)cc1"  # phenol
    analogues = generate_analogues(test_smiles, 10)
    check("Generates analogues", len(analogues) > 0, f"{len(analogues)} generated")
    check("No duplicates", len(analogues) == len(set(analogues)))
    check("Base compound excluded", test_smiles not in analogues)

    # All analogues are valid SMILES
    from rdkit import Chem
    all_valid = all(Chem.MolFromSmiles(s) is not None for s in analogues)
    check("All analogues are valid SMILES", all_valid)

except Exception as e:
    check("Analogue generation import", False, str(e))


# ── 4. Drug-likeness Filter ───────────────────────────────────────────────────
print("\n[4] Drug-likeness Filter (Lipinski + ESOL)")
try:
    from backend.pipeline.drug_likeness import check_drug_likeness, filter_by_drug_likeness

    # Aspirin — should pass
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    result = check_drug_likeness(aspirin)
    check("Aspirin passes Lipinski", result.passed, f"MW={result.mw}")
    check("Aspirin MW correct range", 150 < result.mw < 200, f"{result.mw}")
    check("Solubility class assigned", result.solubility_class in ["High","Medium","Low","Very Low"])

    # Batch filter
    smiles_batch = [aspirin, "C" * 60, "c1ccc(O)cc1"]  # aspirin, huge mol, phenol
    passed, all_results = filter_by_drug_likeness(smiles_batch)
    check("Batch filter runs", len(all_results) == 3)
    check("Filters out oversized molecule", len(passed) < 3)

except Exception as e:
    check("Drug-likeness import", False, str(e))


# ── 5. PDB Fetcher ────────────────────────────────────────────────────────────
print("\n[5] PDB Fetcher (RCSB API)")
try:
    from backend.utils.pdb_fetcher import fetch_pdb_metadata

    # Test metadata fetch (no file download needed)
    async def test_pdb():
        meta = await fetch_pdb_metadata("1HSG")   # HIV protease — classic docking target
        return meta

    meta = asyncio.run(test_pdb())
    check("RCSB metadata fetch", "pdb_id" in meta, meta.get("title","")[:40])
    check("Title returned", bool(meta.get("title")))

except Exception as e:
    check("PDB fetcher", False, str(e))


# ── 6. Binding Site Detection ─────────────────────────────────────────────────
print("\n[6] Binding Site Detection")
try:
    from backend.utils.binding_site import get_binding_site_box
    from backend.models.schemas import BindingSiteMode, BindingSiteCoords

    # Test with explicit coordinates (no PDB file needed)
    center, box_size = get_binding_site_box(
        pdb_path="nonexistent.pdb",  # not used in coord mode
        mode=BindingSiteMode.coordinates,
        coords=BindingSiteCoords(x=10.0, y=20.0, z=5.0, box_size=22.5),
    )
    check("Coordinate mode returns correct center", center == (10.0, 20.0, 5.0))
    check("Coordinate mode returns correct box", box_size == 22.5)

except Exception as e:
    check("Binding site detection", False, str(e))


# ── 7. Scorer ─────────────────────────────────────────────────────────────────
print("\n[7] Scorer & Ranker")
try:
    from backend.pipeline.scorer import score_and_rank, get_score_breakdown
    from backend.models.schemas import (
        CompoundResult, LipinskiResult, ADMETResult,
        DockingResult, DockingPose, BindingPrefilterResult
    )

    # Build two mock CompoundResults with different quality
    c1 = CompoundResult(
        smiles="c1ccc(O)cc1",
        canonical_smiles="Oc1ccccc1",
        lipinski=LipinskiResult(
            passed=True, mw=94.1, logp=1.5, hbd=1, hba=1,
            logs=-1.2, solubility_class="High"
        ),
        admet=ADMETResult(
            passed=True, herg_inhibition=0.2, caco2_permeability=-4.5,
            bbb_penetration=0.4, hepatotoxicity=0.1,
            oral_bioavailability=0.7, flags=[]
        ),
        docking=DockingResult(
            passed=True, best_affinity_kcal=-9.2, cnn_score=0.85,
            poses=[DockingPose(rank=1, affinity_kcal=-9.2, rmsd_lb=0.0, rmsd_ub=0.0)]
        ),
        binding_prefilter=BindingPrefilterResult(
            passed=True, predicted_affinity_kcal=-8.5, confidence=0.78
        ),
    )
    c2 = CompoundResult(
        smiles="CCO",
        canonical_smiles="CCO",
        lipinski=LipinskiResult(
            passed=True, mw=46.0, logp=-0.3, hbd=1, hba=1,
            logs=0.5, solubility_class="High"
        ),
        admet=ADMETResult(
            passed=False, herg_inhibition=0.7, caco2_permeability=-6.0,
            bbb_penetration=0.8, hepatotoxicity=0.6,
            oral_bioavailability=0.2,
            flags=["High hERG cardiotoxicity risk", "Hepatotoxic"]
        ),
        docking=DockingResult(
            passed=False, best_affinity_kcal=-4.1, cnn_score=0.22,
            poses=[DockingPose(rank=1, affinity_kcal=-4.1, rmsd_lb=0.0, rmsd_ub=0.0)]
        ),
    )

    ranked = score_and_rank([c2, c1])  # intentionally give c1 second
    check("Scorer assigns ranks", all(c.rank is not None for c in ranked))
    check("Better compound ranks first", ranked[0].canonical_smiles == "Oc1ccccc1")
    check("Scores are 0–100", all(0 <= c.final_score <= 100 for c in ranked))
    check("Score difference is meaningful", ranked[0].final_score > ranked[1].final_score)

    breakdown = get_score_breakdown(ranked[0])
    check("Score breakdown has docking", "docking_affinity" in breakdown)
    check("Score breakdown has ADMET", "admet_safety" in breakdown)
    check("Final score in breakdown", "final_score" in breakdown)

except Exception as e:
    check("Scorer", False, str(e))


# ── 8. FastAPI app import ─────────────────────────────────────────────────────
print("\n[8] FastAPI App Import")
try:
    from backend.main import app
    from fastapi.testclient import TestClient
    client = TestClient(app)

    resp = client.get("/")
    check("Root endpoint responds", resp.status_code == 200)
    check("Returns version", resp.json().get("version") == "1.0.0")

    resp2 = client.get("/docs")
    check("Swagger UI accessible", resp2.status_code == 200)

except Exception as e:
    check("FastAPI app import", False, str(e))


# ── 9. ADMET-AI (optional — heavy model) ─────────────────────────────────────
print("\n[9] ADMET-AI (may be slow on first run — downloads model weights)")
try:
    from admet_ai import ADMETModel
    model = ADMETModel()
    df = model.predict(smiles=["c1ccccc1", "CC(=O)O"])
    check("ADMET-AI predicts batch", len(df) == 2)
    check("hERG column exists", "hERG" in df.columns or any("hERG" in c for c in df.columns))
except ImportError:
    print(f"  {WARN} ADMET-AI not installed — run: pip install admet-ai")
    results.append(None)
except Exception as e:
    check("ADMET-AI inference", False, str(e))


# ── Summary ───────────────────────────────────────────────────────────────────
print("\n" + "="*55)
hard_results = [r for r in results if r is not None]
passed_count = sum(hard_results)
total_count  = len(hard_results)
all_passed   = all(hard_results)

print(f"  Result: {passed_count}/{total_count} checks passed")
if all_passed:
    print(f"  {PASS} ALL CHECKS PASSED — backend logic is solid")
else:
    failed = total_count - passed_count
    print(f"  {FAIL} {failed} check(s) failed — see details above")
print("="*55 + "\n")
