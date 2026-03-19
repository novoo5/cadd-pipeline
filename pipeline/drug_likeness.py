from typing import List, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from backend.models.schemas import LipinskiResult
from backend.utils.smiles_validator import smiles_to_mol


# ── Lipinski Rule of 5 hard limits ──────────────────────────────────────────
LIPINSKI_MW_MAX       = 500.0   # Molecular weight ≤ 500 Da
LIPINSKI_LOGP_MAX     = 5.0     # LogP ≤ 5
LIPINSKI_HBD_MAX      = 5       # H-bond donors ≤ 5
LIPINSKI_HBA_MAX      = 10      # H-bond acceptors ≤ 10

# One violation is allowed (Lipinski's original rule)
LIPINSKI_MAX_VIOLATIONS = 1

# ESOL solubility thresholds (LogS, mol/L)
SOLUBILITY_HIGH      =  0.0     # LogS ≥ 0    → High
SOLUBILITY_MEDIUM    = -2.0     # LogS ≥ -2   → Medium
SOLUBILITY_LOW       = -4.0     # LogS ≥ -4   → Low
                                 # LogS <  -4  → Very Low


def check_drug_likeness(smiles: str) -> LipinskiResult:
    """
    Runs Lipinski Rule of 5 + ESOL solubility prediction on a SMILES string.
    Returns a fully populated LipinskiResult.
    """
    mol = smiles_to_mol(smiles)

    mw    = Descriptors.MolWt(mol)
    logp  = Descriptors.MolLogP(mol)
    hbd   = rdMolDescriptors.CalcNumHBD(mol)
    hba   = rdMolDescriptors.CalcNumHBA(mol)
    logs  = _calc_esol_logs(mol)
    sol_class = _solubility_class(logs)

    violations = sum([
        mw   > LIPINSKI_MW_MAX,
        logp > LIPINSKI_LOGP_MAX,
        hbd  > LIPINSKI_HBD_MAX,
        hba  > LIPINSKI_HBA_MAX,
    ])
    passed = violations <= LIPINSKI_MAX_VIOLATIONS

    return LipinskiResult(
        passed=passed,
        mw=round(mw, 3),
        logp=round(logp, 3),
        hbd=hbd,
        hba=hba,
        logs=round(logs, 3),
        solubility_class=sol_class,
    )


def filter_by_drug_likeness(smiles_list: List[str]) -> Tuple[List[str], List[LipinskiResult]]:
    """
    Batch filter: runs check_drug_likeness on every SMILES in the list.
    Returns:
      - passed_smiles: only compounds where LipinskiResult.passed == True
      - all_results: LipinskiResult for every input compound (for reporting)
    """
    passed_smiles: List[str] = []
    all_results: List[LipinskiResult] = []

    for smi in smiles_list:
        try:
            result = check_drug_likeness(smi)
            all_results.append(result)
            if result.passed:
                passed_smiles.append(smi)
        except Exception as e:
            # If RDKit can't handle it, skip it
            all_results.append(LipinskiResult(
                passed=False,
                mw=0.0, logp=0.0, hbd=0, hba=0,
                logs=-999.0,
                solubility_class="Invalid",
            ))

    return passed_smiles, all_results


def _calc_esol_logs(mol) -> float:
    """
    ESOL (Delaney) model for aqueous solubility estimation.
    Formula: LogS = 0.16 - 0.63*clogP - 0.0062*MW + 0.066*RB - 0.74*AP
    where:
      clogP = calculated LogP
      MW    = molecular weight
      RB    = rotatable bond count
      AP    = aromatic proportion (aromatic atoms / heavy atoms)

    Reference: Delaney, J.S. (2004) JCIM 44(3):1000-1005
    """
    mw    = Descriptors.MolWt(mol)
    clogp = Descriptors.MolLogP(mol)
    rb    = rdMolDescriptors.CalcNumRotatableBonds(mol)

    heavy_atoms    = mol.GetNumHeavyAtoms()
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    ap = aromatic_atoms / heavy_atoms if heavy_atoms > 0 else 0.0

    logs = 0.16 - (0.63 * clogp) - (0.0062 * mw) + (0.066 * rb) - (0.74 * ap)
    return logs


def _solubility_class(logs: float) -> str:
    if logs >= SOLUBILITY_HIGH:
        return "High"
    elif logs >= SOLUBILITY_MEDIUM:
        return "Medium"
    elif logs >= SOLUBILITY_LOW:
        return "Low"
    else:
        return "Very Low"
