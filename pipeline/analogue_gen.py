import random
from typing import List, Set
from rdkit import Chem
from rdkit.Chem import BRICS, AllChem, Descriptors
from backend.utils.smiles_validator import smiles_to_mol, mol_to_smiles


def generate_analogues(base_smiles: str, num_analogues: int) -> List[str]:
    """
    Generates structurally diverse analogues of the base compound using three strategies:
      1. BRICS fragmentation + recombination (primary — most chemically valid)
      2. Single-atom mutations (swap common functional groups)
      3. Ring-aware substituent addition

    Returns a deduplicated list of canonical SMILES strings (excluding the base compound).
    Always attempts to return exactly num_analogues, but may return fewer if
    chemical space around the scaffold is limited.
    """
    base_mol = smiles_to_mol(base_smiles)
    base_canonical = mol_to_smiles(base_mol)

    collected: Set[str] = set()
    collected.add(base_canonical)  # exclude base from results

    # --- Strategy 1: BRICS fragmentation + recombination ---
    brics_results = _brics_recombination(base_mol, num_analogues * 2)
    for smi in brics_results:
        if len(collected) >= num_analogues + 1:
            break
        collected.add(smi)

    # --- Strategy 2: Functional group mutations ---
    if len(collected) - 1 < num_analogues:
        mutation_results = _functional_group_mutations(base_mol, num_analogues)
        for smi in mutation_results:
            if len(collected) >= num_analogues + 1:
                break
            collected.add(smi)

    # --- Strategy 3: Substituent additions ---
    if len(collected) - 1 < num_analogues:
        addition_results = _substituent_additions(base_mol, num_analogues)
        for smi in addition_results:
            if len(collected) >= num_analogues + 1:
                break
            collected.add(smi)

    # Remove base compound from results before returning
    collected.discard(base_canonical)
    return list(collected)[:num_analogues]


def _brics_recombination(mol, target_count: int) -> List[str]:
    """
    Breaks the molecule at BRICS bonds, then recombines fragments randomly.
    BRICS ensures breaks happen at chemically meaningful bond types.
    """
    results = []
    try:
        fragments = list(BRICS.BRICSDecompose(mol))
        if len(fragments) < 2:
            return results

        # Build fragment mols, skip invalid ones
        frag_mols = []
        for f in fragments:
            m = Chem.MolFromSmiles(f)
            if m is not None:
                frag_mols.append(m)

        if len(frag_mols) < 2:
            return results

        # Generate recombined molecules
        attempts = 0
        max_attempts = target_count * 10
        while len(results) < target_count and attempts < max_attempts:
            attempts += 1
            try:
                # Randomly pick 2 fragments and recombine
                chosen = random.sample(frag_mols, min(2, len(frag_mols)))
                new_mols = list(BRICS.BRICSBuild(chosen))
                for new_mol in new_mols:
                    if new_mol is None:
                        continue
                    smi = mol_to_smiles(new_mol)
                    if _is_valid_analogue(new_mol):
                        results.append(smi)
                    if len(results) >= target_count:
                        break
            except Exception:
                continue

    except Exception:
        pass

    return results


# Functional group swap pairs: (SMARTS_to_find, replacement_SMILES_fragment)
_FG_MUTATIONS = [
    ("[OH]",       "OC"),          # hydroxyl → methoxy
    ("[OH]",       "N"),           # hydroxyl → amine
    ("[OH]",       "F"),           # hydroxyl → fluorine (bioisostere)
    ("[NH2]",      "O"),           # amine → hydroxyl
    ("[NH2]",      "NC"),          # amine → methylamine
    ("[F]",        "Cl"),          # F → Cl
    ("[F]",        "Br"),          # F → Br
    ("[Cl]",       "F"),           # Cl → F
    ("[CH3]",      "CC"),          # methyl → ethyl
    ("[CH3]",      "CF"),          # methyl → fluoromethyl
    ("C(=O)O",     "C(=O)N"),      # ester → amide
    ("C(=O)N",     "C(=O)O"),      # amide → ester
]


def _functional_group_mutations(mol, target_count: int) -> List[str]:
    """
    Applies pairwise functional group swaps using SMARTS patterns.
    Each successful swap produces one new analogue.
    """
    results = []
    base_smi = mol_to_smiles(mol)

    for smarts, replacement in _FG_MUTATIONS:
        if len(results) >= target_count:
            break
        try:
            pattern = Chem.MolFromSmarts(smarts)
            repl_mol = Chem.MolFromSmiles(replacement)
            if pattern is None or repl_mol is None:
                continue

            # AllChem.ReplaceSubstructs returns a tuple of product mols
            products = AllChem.ReplaceSubstructs(mol, pattern, repl_mol)
            for product in products:
                if product is None:
                    continue
                try:
                    smi = mol_to_smiles(product)
                    if smi != base_smi and _is_valid_analogue(product):
                        results.append(smi)
                except Exception:
                    continue

        except Exception:
            continue

    return results


# Substituents to try adding at open positions
_SUBSTITUENTS = [
    "F", "Cl", "Br",
    "C", "CC", "CCC",
    "OC", "OCC",
    "N", "NC", "NCC",
    "C(F)(F)F",      # trifluoromethyl
    "C#N",           # nitrile
    "S(=O)(=O)N",    # sulfonamide
    "C(=O)N",        # amide
]


def _substituent_additions(mol, target_count: int) -> List[str]:
    """
    Attempts to add substituents to aromatic carbons with open hydrogen positions.
    Uses atom mapping to find substitutable positions.
    """
    results = []
    base_smi = mol_to_smiles(mol)

    rw_mol = Chem.RWMol(mol)
    aromatic_c_with_h = [
        atom.GetIdx()
        for atom in rw_mol.GetAtoms()
        if atom.GetAtomicNum() == 6
        and atom.GetIsAromatic()
        and atom.GetTotalNumHs() > 0
    ]

    for idx in aromatic_c_with_h:
        if len(results) >= target_count:
            break
        for sub_smi in _SUBSTITUENTS:
            if len(results) >= target_count:
                break
            try:
                sub_mol = Chem.MolFromSmiles(sub_smi)
                if sub_mol is None:
                    continue
                pattern = Chem.MolFromSmarts("[cH]")
                if pattern is None:
                    continue
                products = AllChem.ReplaceSubstructs(mol, pattern, sub_mol, replaceAll=False)
                for product in products:
                    if product is None:
                        continue
                    try:
                        smi = mol_to_smiles(product)
                        if smi != base_smi and _is_valid_analogue(product):
                            results.append(smi)
                    except Exception:
                        continue
            except Exception:
                continue

    return results


def _is_valid_analogue(mol) -> bool:
    """
    Fast sanity check before accepting an analogue.
    Rejects molecules that are too small, too large, or have impossible valence.
    """
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
        num_atoms = mol.GetNumAtoms()
        if num_atoms < 5 or num_atoms > 200:
            return False
        mw = Descriptors.MolWt(mol)
        if mw < 100 or mw > 1500:
            return False
        return True
    except Exception:
        return False
