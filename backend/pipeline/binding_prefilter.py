from typing import List, Tuple
from backend.models.schemas import BindingPrefilterResult


try:
    import deepchem as dc
    import numpy as np
    _DC_AVAILABLE = True
except ImportError:
    _DC_AVAILABLE = False


_MODEL = None


def _load_model():
    global _MODEL
    if _MODEL is not None:
        return _MODEL
    if not _DC_AVAILABLE:
        return None
    try:
        _MODEL = dc.models.AttentiveFPModel(
            mode="regression",
            n_tasks=1,
            model_dir=None,
        )
        return _MODEL
    except Exception as e:
        print(f"[Prefilter] Could not load AttentiveFP model: {e}")
        return None


def run_binding_prefilter(
    smiles_list: List[str],
    top_n: int = 10,
) -> Tuple[List[str], List[BindingPrefilterResult]]:
    if not smiles_list:
        return [], []

    if not _DC_AVAILABLE:
        return _fallback_prefilter(smiles_list, top_n)

    model = _load_model()
    if model is None:
        return _fallback_prefilter(smiles_list, top_n)

    try:
        featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
        features = featurizer.featurize(smiles_list)
        dataset = dc.data.NumpyDataset(X=features)
        predictions = model.predict(dataset).flatten()

        results = []
        scored = []

        for i, smi in enumerate(smiles_list):
            pred = float(predictions[i]) if i < len(predictions) else 0.0
            affinity_approx = -1.0 * abs(pred)
            confidence = min(abs(pred) / 10.0, 1.0)

            result = BindingPrefilterResult(
                passed=False,
                predicted_affinity_kcal=round(affinity_approx, 4),
                confidence=round(confidence, 4),
            )
            results.append(result)
            scored.append((i, affinity_approx))

        # More negative = better predicted binder
        scored.sort(key=lambda x: x[1])
        actual_top_n = min(top_n, len(smiles_list))
        top_indices = {idx for idx, _ in scored[:actual_top_n]}

        top_smiles = []
        for i, result in enumerate(results):
            if i in top_indices:
                result.passed = True
                top_smiles.append(smiles_list[i])

        return top_smiles, results

    except Exception as e:
        print(f"[Prefilter] DeepChem inference failed: {e}. Using fallback.")
        return _fallback_prefilter(smiles_list, top_n)


def _fallback_prefilter(
    smiles_list: List[str],
    top_n: int,
) -> Tuple[List[str], List[BindingPrefilterResult]]:
    """
    When DeepChem unavailable, pass ALL compounds through only up to top_n.
    Docking is the real filter — prefilter is just a compute budget gate.
    """
    actual_top = min(top_n, len(smiles_list))
    results = []

    for i, smi in enumerate(smiles_list):
        passed = i < actual_top
        results.append(
            BindingPrefilterResult(
                passed=passed,
                predicted_affinity_kcal=-6.0,
                confidence=0.0,
            )
        )

    top_smiles = smiles_list[:actual_top]
    return top_smiles, results