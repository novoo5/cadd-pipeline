import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "backend"))

from pipeline.protein_prep import prepare_ligand

print("Trying prep...")
try:
    path = prepare_ligand('CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'test_job', 0)
    print(f"Prep success! Saved to {path}")
except Exception as e:
    import traceback
    print("CRASHED:")
    print(traceback.format_exc())
