"""
Exhaustive analogue generator supporting up to 10000 compounds.

Strategy stack (applied in order, then looped with 2nd-order mutations):
  1. BRICS fragmentation + recombination  (chemically valid, diverse)
  2. Functional group swaps               (bioisostere replacements)
  3. Substituent additions                (halogen scan, OH, CN, CF3, OMe, NH2)
  4. Direct atom mutations                (C/N/O/S swaps, aza-scan on rings)
  5. Side-chain edits                     (methyl add/remove, chain extend/shorten)

Solubility filter (Delaney ESOL) is applied at generation time when requested.
modes: "soluble_only" (logS > -3) | "allow_slightly" (logS > -5) | "all" (no filter)
"""

import random
from typing import List, Set, Optional, Callable  # ← added Callable
from rdkit import Chem
from rdkit.Chem import BRICS, AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.rdchem import RWMol, BondType
from backend.utils.smiles_validator import smiles_to_mol, mol_to_smiles


# ── Solubility (Delaney ESOL) ─────────────────────────────────────────────────

def esol_logS(mol) -> float:
    """
    Delaney ESOL model:
      logS = 0.16 − 0.63·cLogP − 0.0062·MW + 0.066·nHBA − 0.74·nRB
    More negative = less soluble.
    """
    try:
        logP = Descriptors.MolLogP(mol)
        mw   = Descriptors.MolWt(mol)
        hba  = rdMolDescriptors.CalcNumHBA(mol)
        rb   = rdMolDescriptors.CalcNumRotatableBonds(mol)
        return 0.16 - 0.63 * logP - 0.0062 * mw + 0.066 * hba - 0.74 * rb
    except Exception:
        return -99.0


def solubility_class(logs: float) -> str:
    if logs > -1:  return "highly_soluble"
    if logs > -3:  return "soluble"
    if logs > -5:  return "slightly_soluble"
    return "insoluble"


_SOLUBILITY_THRESHOLDS = {
    "soluble_only":   -3.0,
    "allow_slightly": -5.0,
    "all":            -99.0,
}


# ── Validity gate ─────────────────────────────────────────────────────────────

def _is_valid(
    mol,
    base_canonical: str,
    seen:           Set[str],
    logs_threshold: float,
    min_mw:         float,
    max_mw:         float,
) -> Optional[str]:
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        smi = Chem.MolToSmiles(mol)
        if not smi or smi in seen or smi == base_canonical:
            return None
        mw = Descriptors.MolWt(mol)
        if not (min_mw <= mw <= max_mw):
            return None
        na = mol.GetNumAtoms()
        if na < 5 or na > 200:
            return None
        if logs_threshold > -90.0 and esol_logS(mol) < logs_threshold:
            return None
        return smi
    except Exception:
        return None


# ── Strategy 1: BRICS ─────────────────────────────────────────────────────────

def _brics_batch(mol, target: int, base_canonical: str, seen: Set[str],
                 logs_threshold: float, min_mw: float, max_mw: float) -> List[str]:
    results = []
    try:
        fragments = list(BRICS.BRICSDecompose(mol))
        if len(fragments) < 2:
            return results
        frag_mols = [Chem.MolFromSmiles(f) for f in fragments]
        frag_mols = [m for m in frag_mols if m is not None]
        if len(frag_mols) < 2:
            return results

        max_attempts = target * 15
        for _ in range(max_attempts):
            if len(results) >= target:
                break
            try:
                chosen   = random.sample(frag_mols, min(2, len(frag_mols)))
                new_mols = list(BRICS.BRICSBuild(chosen))
                for new_mol in new_mols:
                    smi = _is_valid(new_mol, base_canonical, seen, logs_threshold, min_mw, max_mw)
                    if smi:
                        seen.add(smi)
                        results.append(smi)
                        if len(results) >= target:
                            break
            except Exception:
                continue
    except Exception:
        pass
    return results


# ── Strategy 2: Functional group swaps ───────────────────────────────────────

_FG_MUTATIONS = [
    ("[OH]",   "OC"),     ("[OH]",  "N"),      ("[OH]",  "F"),
    ("[NH2]",  "O"),      ("[NH2]", "NC"),      ("[NH2]", "NCC"),
    ("[F]",    "Cl"),     ("[F]",   "Br"),      ("[F]",   "C(F)(F)F"),
    ("[Cl]",   "F"),      ("[Cl]",  "Br"),
    ("[Br]",   "F"),      ("[Br]",  "Cl"),
    ("[CH3]",  "CC"),     ("[CH3]", "CCC"),     ("[CH3]", "CF"),
    ("[CH3]",  "C(F)(F)F"),
    ("C(=O)O", "C(=O)N"), ("C(=O)N", "C(=O)O"),
    ("c[OH]",  "cOC"),    ("c[OH]",  "cN"),
]


def _fg_mutation_batch(mol, target: int, base_canonical: str, seen: Set[str],
                       logs_threshold: float, min_mw: float, max_mw: float) -> List[str]:
    results = []
    mutations = random.sample(_FG_MUTATIONS, len(_FG_MUTATIONS))
    for smarts, replacement in mutations:
        if len(results) >= target:
            break
        try:
            pattern  = Chem.MolFromSmarts(smarts)
            repl_mol = Chem.MolFromSmiles(replacement)
            if pattern is None or repl_mol is None:
                continue
            products = AllChem.ReplaceSubstructs(mol, pattern, repl_mol)
            for product in products:
                smi = _is_valid(product, base_canonical, seen, logs_threshold, min_mw, max_mw)
                if smi:
                    seen.add(smi)
                    results.append(smi)
        except Exception:
            continue
    return results


# ── Strategy 3: Substituent additions ────────────────────────────────────────

_SUBSTITUENTS = [
    "F", "Cl", "Br",
    "C", "CC", "CCC",
    "OC", "OCC", "OCCC",
    "N", "NC", "NCC",
    "C(F)(F)F",
    "C#N",
    "S(=O)(=O)N",
    "C(=O)N",
    "C(=O)O",
    "OC(F)(F)F",
    "SC",
]


def _substituent_batch(mol, target: int, base_canonical: str, seen: Set[str],
                       logs_threshold: float, min_mw: float, max_mw: float) -> List[str]:
    results = []
    pattern = Chem.MolFromSmarts("[cH]")
    if pattern is None:
        return results
    subs = random.sample(_SUBSTITUENTS, len(_SUBSTITUENTS))
    for sub_smi in subs:
        if len(results) >= target:
            break
        try:
            sub_mol = Chem.MolFromSmiles(sub_smi)
            if sub_mol is None:
                continue
            products = AllChem.ReplaceSubstructs(mol, pattern, sub_mol, replaceAll=False)
            for product in products:
                smi = _is_valid(product, base_canonical, seen, logs_threshold, min_mw, max_mw)
                if smi:
                    seen.add(smi)
                    results.append(smi)
                    if len(results) >= target:
                        break
        except Exception:
            continue
    return results


# ── Strategy 4: Direct atom mutations (RWMol) ─────────────────────────────────

_ATOM_SUBS = {
    "C": ["N", "O", "S"],
    "N": ["C", "O"],
    "O": ["N", "S"],
    "S": ["O", "N"],
}
_PERIODIC = Chem.GetPeriodicTable()


def _rwmol_mutations_batch(mol, target: int, base_canonical: str, seen: Set[str],
                           logs_threshold: float, min_mw: float, max_mw: float) -> List[str]:
    results = []
    atoms = list(mol.GetAtoms())

    def _try_rw(func):
        rw = RWMol(mol)
        try:
            func(rw)
            Chem.SanitizeMol(rw)
            smi = _is_valid(rw.GetMol(), base_canonical, seen, logs_threshold, min_mw, max_mw)
            if smi:
                seen.add(smi)
                results.append(smi)
        except Exception:
            pass

    random.shuffle(atoms)
    for atom in atoms:
        if len(results) >= target:
            break
        sym = atom.GetSymbol()
        idx = atom.GetIdx()

        if sym in _ATOM_SUBS and not atom.IsInRing():
            rep = random.choice(_ATOM_SUBS[sym])
            def _sub(rw, _idx=idx, _rep=rep):
                rw.GetAtomWithIdx(_idx).SetAtomicNum(_PERIODIC.GetAtomicNumber(_rep))
            _try_rw(_sub)

        if sym == "C" and atom.GetIsAromatic() and atom.IsInRing():
            def _aza(rw, _idx=idx):
                rw.GetAtomWithIdx(_idx).SetAtomicNum(7)
            _try_rw(_aza)

        if sym in ("F", "Cl", "Br"):
            hals = [x for x in ("F", "Cl", "Br") if x != sym]
            rep  = random.choice(hals)
            def _hal(rw, _idx=idx, _rep=rep):
                rw.GetAtomWithIdx(_idx).SetAtomicNum(_PERIODIC.GetAtomicNumber(_rep))
            _try_rw(_hal)

        if sym == "C" and atom.GetIsAromatic() and atom.GetTotalNumHs() > 0:
            hal = random.choice(["F", "Cl", "Br"])
            def _hadd(rw, _idx=idx, _hal=hal):
                new = rw.AddAtom(Chem.Atom(_hal))
                rw.AddBond(_idx, new, BondType.SINGLE)
            _try_rw(_hadd)

    return results


# ── Strategy 5: Side-chain edits ──────────────────────────────────────────────

def _chain_edits_batch(mol, target: int, base_canonical: str, seen: Set[str],
                       logs_threshold: float, min_mw: float, max_mw: float) -> List[str]:
    results = []
    atoms = list(mol.GetAtoms())
    random.shuffle(atoms)

    def _try_rw(func):
        rw = RWMol(mol)
        try:
            func(rw)
            Chem.SanitizeMol(rw)
            smi = _is_valid(rw.GetMol(), base_canonical, seen, logs_threshold, min_mw, max_mw)
            if smi:
                seen.add(smi)
                results.append(smi)
        except Exception:
            pass

    for atom in atoms:
        if len(results) >= target:
            break
        sym = atom.GetSymbol()
        idx = atom.GetIdx()

        if sym in ("C", "N") and not atom.GetIsAromatic() and atom.GetTotalNumHs() > 0:
            def _madd(rw, _idx=idx):
                c = rw.AddAtom(Chem.Atom("C"))
                rw.AddBond(_idx, c, BondType.SINGLE)
            _try_rw(_madd)

        if sym == "C" and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 3:
            def _mrem(rw, _idx=idx):
                rw.RemoveAtom(_idx)
            _try_rw(_mrem)

        if sym == "C" and not atom.GetIsAromatic() and atom.GetDegree() <= 3 and atom.GetTotalNumHs() > 0:
            def _cext(rw, _idx=idx):
                c = rw.AddAtom(Chem.Atom("C"))
                rw.AddBond(_idx, c, BondType.SINGLE)
            _try_rw(_cext)

        if sym == "C" and not atom.GetIsAromatic() and atom.GetDegree() == 1 and atom.GetTotalNumHs() >= 2:
            def _csh(rw, _idx=idx):
                rw.RemoveAtom(_idx)
            _try_rw(_csh)

        if sym == "C" and atom.GetIsAromatic() and atom.GetTotalNumHs() > 0:
            def _ohadd(rw, _idx=idx):
                o = rw.AddAtom(Chem.Atom("O"))
                rw.AddBond(_idx, o, BondType.SINGLE)
            _try_rw(_ohadd)

    return results


# ── Public API ────────────────────────────────────────────────────────────────

def generate_analogues(
    base_smiles:             str,
    num_analogues:           int,
    filter_solubility:       str                              = "all",
    min_mw:                  float                           = 100.0,
    max_mw:                  float                           = 800.0,
    max_attempts_multiplier: int                             = 40,
    on_progress:             Optional[Callable[[int, int], None]] = None,  # ← NEW
) -> List[str]:
    """
    Generate up to num_analogues unique analogues of base_smiles.

    Args:
        base_smiles:             Input SMILES.
        num_analogues:           Target count (1–10000).
        filter_solubility:       "soluble_only" | "allow_slightly" | "all"
        min_mw / max_mw:         MW window for acceptance.
        max_attempts_multiplier: Total strategy loop budget = n × multiplier.
        on_progress:             Optional callback(current, total) fired every
                                 10 new compounds — used for real-time UI updates.

    Returns:
        List of canonical SMILES (up to num_analogues, may be fewer if
        chemical space is genuinely exhausted).
    """
    base_mol = smiles_to_mol(base_smiles)
    if base_mol is None:
        return []

    base_canonical     = mol_to_smiles(base_mol)
    logs_threshold     = _SOLUBILITY_THRESHOLDS.get(filter_solubility, -99.0)
    seen: Set[str]     = {base_canonical}
    results: List[str] = []

    # ── Scale multiplier down at high counts to keep runtime sane ────────────
    if num_analogues > 5000:
        max_attempts_multiplier = 8
    elif num_analogues > 2000:
        max_attempts_multiplier = 12
    elif num_analogues > 1000:
        max_attempts_multiplier = 20
    elif num_analogues > 500:
        max_attempts_multiplier = 25
    elif num_analogues > 200:
        max_attempts_multiplier = 30
    # ≤200 keeps the default 40×

    # ── Progress reporter: fires every 10 new compounds ───────────────────────
    _last_reported = [0]

    def _maybe_report():
        if on_progress is None:
            return
        if len(results) - _last_reported[0] >= 10 or len(results) == num_analogues:
            on_progress(len(results), num_analogues)
            _last_reported[0] = len(results)

    # ── Pass 1: all 5 strategies on the base molecule ─────────────────────────
    batch = _brics_batch(base_mol, num_analogues - len(results), base_canonical, seen, logs_threshold, min_mw, max_mw)
    results.extend(batch)
    _maybe_report()

    if len(results) < num_analogues:
        batch = _fg_mutation_batch(base_mol, num_analogues - len(results), base_canonical, seen, logs_threshold, min_mw, max_mw)
        results.extend(batch)
        _maybe_report()

    if len(results) < num_analogues:
        batch = _substituent_batch(base_mol, num_analogues - len(results), base_canonical, seen, logs_threshold, min_mw, max_mw)
        results.extend(batch)
        _maybe_report()

    if len(results) < num_analogues:
        batch = _rwmol_mutations_batch(base_mol, num_analogues - len(results), base_canonical, seen, logs_threshold, min_mw, max_mw)
        results.extend(batch)
        _maybe_report()

    if len(results) < num_analogues:
        batch = _chain_edits_batch(base_mol, num_analogues - len(results), base_canonical, seen, logs_threshold, min_mw, max_mw)
        results.extend(batch)
        _maybe_report()

    # ── Pass 2+: 2nd-order mutations from accepted analogues ──────────────────
    budget = num_analogues * max_attempts_multiplier
    iters  = 0

    while len(results) < num_analogues and iters < budget:
        iters += 1
        if not results:
            break

        parent_smi = random.choice(results)
        parent_mol = Chem.MolFromSmiles(parent_smi)
        if parent_mol is None:
            continue

        strategy = random.randint(1, 5)
        need     = num_analogues - len(results)

        if strategy == 1:
            batch = _brics_batch(parent_mol, min(need, 5), base_canonical, seen, logs_threshold, min_mw, max_mw)
        elif strategy == 2:
            batch = _fg_mutation_batch(parent_mol, min(need, 5), base_canonical, seen, logs_threshold, min_mw, max_mw)
        elif strategy == 3:
            batch = _substituent_batch(parent_mol, min(need, 5), base_canonical, seen, logs_threshold, min_mw, max_mw)
        elif strategy == 4:
            batch = _rwmol_mutations_batch(parent_mol, min(need, 5), base_canonical, seen, logs_threshold, min_mw, max_mw)
        else:
            batch = _chain_edits_batch(parent_mol, min(need, 5), base_canonical, seen, logs_threshold, min_mw, max_mw)

        results.extend(batch)
        _maybe_report()

    return results[:num_analogues]


def get_solubility_info(smiles: str) -> dict:
    """Return ESOL logS and solubility class for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"logs": None, "class": "unknown"}
    logs = esol_logS(mol)
    return {"logs": round(logs, 3), "class": solubility_class(logs)}
