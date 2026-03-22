from typing import List, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from backend.models.schemas import LipinskiResult
from backend.utils.smiles_validator import smiles_to_mol


# Lipinski defaults (overridable via request params)
DEFAULT_MW_MIN              = 200.0
DEFAULT_MW_MAX              = 500.0
LIPINSKI_LOGP_MAX           = 4.5
LIPINSKI_HBD_MAX            = 5
LIPINSKI_HBA_MAX            = 10
DEFAULT_MAX_VIOLATIONS      = 1   # classic strict Lipinski

# ESOL solubility thresholds (LogS, mol/L)
SOLUBILITY_HIGH   =  0.0
SOLUBILITY_MEDIUM = -2.0
SOLUBILITY_LOW    = -4.0


def check_drug_likeness(
    smiles: str,
    mw_min: float = DEFAULT_MW_MIN,
    mw_max: float = DEFAULT_MW_MAX,
    max_violations: Optional[int] = DEFAULT_MAX_VIOLATIONS,
) -> LipinskiResult:
    """
    Evaluates Lipinski RO5 + ESOL for a single SMILES.
    max_violations=None → "Ignore completely": passed is always True,
    violation counting still runs so the values are stored for display.
    """
    mol = smiles_to_mol(smiles)

    mw        = Descriptors.MolWt(mol)
    logp      = Descriptors.MolLogP(mol)
    hbd       = rdMolDescriptors.CalcNumHBD(mol)
    hba       = rdMolDescriptors.CalcNumHBA(mol)
    logs      = _calc_esol_logs(mol)
    sol_class = _solubility_class(logs)

    violations = sum([
        mw   < mw_min,
        mw   > mw_max,
        logp > LIPINSKI_LOGP_MAX,
        hbd  > LIPINSKI_HBD_MAX,
        hba  > LIPINSKI_HBA_MAX,
    ])

    # None = "Ignore completely" → always pass through
    if max_violations is None:
        passed = True
    else:
        passed = violations <= max_violations

    return LipinskiResult(
        passed=passed,
        mw=round(mw, 3),
        logp=round(logp, 3),
        hbd=hbd,
        hba=hba,
        logs=round(logs, 3),
        solubility_class=sol_class,
    )


def filter_by_drug_likeness(
    smiles_list: List[str],
    mw_min: float = DEFAULT_MW_MIN,
    mw_max: float = DEFAULT_MW_MAX,
    max_violations: Optional[int] = DEFAULT_MAX_VIOLATIONS,
) -> Tuple[List[str], List[LipinskiResult]]:
    """
    Filters a list of SMILES by drug-likeness.
    max_violations=None → skip filter, all valid compounds pass through.
    """
    passed_smiles: List[str] = []
    all_results: List[LipinskiResult] = []

    for smi in smiles_list:
        try:
            result = check_drug_likeness(
                smi,
                mw_min=mw_min,
                mw_max=mw_max,
                max_violations=max_violations,
            )
            all_results.append(result)
            if result.passed:
                passed_smiles.append(smi)
        except Exception:
            all_results.append(LipinskiResult(
                passed=False,
                mw=0.0, logp=0.0, hbd=0, hba=0,
                logs=-999.0,
                solubility_class="Invalid",
            ))

    return passed_smiles, all_results


def _calc_esol_logs(mol) -> float:
    """
    ESOL (Delaney) model for aqueous solubility.
    LogS = 0.16 - 0.63*clogP - 0.0062*MW + 0.066*RB - 0.74*AP
    Reference: Delaney, J.S. (2004) JCIM 44(3):1000-1005
    """
    mw    = Descriptors.MolWt(mol)
    clogp = Descriptors.MolLogP(mol)
    rb    = rdMolDescriptors.CalcNumRotatableBonds(mol)

    heavy_atoms    = mol.GetNumHeavyAtoms()
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    ap = aromatic_atoms / heavy_atoms if heavy_atoms > 0 else 0.0

    return 0.16 - (0.63 * clogp) - (0.0062 * mw) + (0.066 * rb) - (0.74 * ap)


def _solubility_class(logs: float) -> str:
    if logs >= SOLUBILITY_HIGH:
        return "High"
    elif logs >= SOLUBILITY_MEDIUM:
        return "Medium"
    elif logs >= SOLUBILITY_LOW:
        return "Low"
    else:
        return "Very Low"
