from rdkit import Chem
from rdkit.Chem import Descriptors
from typing import Tuple, Optional


def validate_and_canonicalize(smiles: str) -> Tuple[bool, Optional[str], Optional[str]]:
    """
    Validates a SMILES string using RDKit.
    Returns: (is_valid, canonical_smiles, error_message)
    """
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return False, None, "RDKit could not parse this SMILES string. Check for typos."

        canonical = Chem.MolToSmiles(mol, canonical=True)

        # Basic sanity checks
        num_atoms = mol.GetNumAtoms()
        if num_atoms < 3:
            return False, None, f"Molecule too small ({num_atoms} atoms). Minimum is 3."
        if num_atoms > 200:
            return False, None, f"Molecule too large ({num_atoms} atoms). Max supported is 200 atoms."

        mw = Descriptors.MolWt(mol)
        if mw > 2000:
            return False, None, f"Molecular weight ({mw:.1f} Da) exceeds 2000 Da hard limit."

        return True, canonical, None

    except Exception as e:
        return False, None, f"Unexpected error during SMILES validation: {str(e)}"


def smiles_to_mol(smiles: str):
    """
    Returns an RDKit Mol object or raises ValueError.
    Always canonicalizes before returning.
    """
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError(f"Invalid SMILES: '{smiles}'")
    return mol


def mol_to_smiles(mol) -> str:
    """Converts RDKit Mol object back to canonical SMILES."""
    return Chem.MolToSmiles(mol, canonical=True)


def get_basic_properties(smiles: str) -> dict:
    """
    Returns a dict of basic molecular properties for a valid SMILES.
    Useful for quick sanity displays before running the full pipeline.
    """
    mol = smiles_to_mol(smiles)
    return {
        "canonical_smiles": mol_to_smiles(mol),
        "num_atoms": mol.GetNumAtoms(),
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "molecular_weight": round(Descriptors.MolWt(mol), 3),
        "formula": _get_formula(mol),
    }


def _get_formula(mol) -> str:
    from rdkit.Chem import rdMolDescriptors
    return rdMolDescriptors.CalcMolFormula(mol)
