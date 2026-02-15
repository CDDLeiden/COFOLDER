"""Minimal SuCOS implementation vendored for COFOLDER.

Adapted conceptually from susanhleung/SuCOS:
- Shape term from protrude distance
- Feature term from RDKit pharmacophore features
- Composite SuCOS = 0.5 * (shape + feature)

License notice:
This file contains adapted logic from susanhleung/SuCOS (MIT License).
See `THIRD_PARTY_LICENSES.md` in the repository root for attribution and
the full MIT license text.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps

KEEP_FEATURES = (
    "Donor",
    "Acceptor",
    "NegIonizable",
    "PosIonizable",
    "ZnBinder",
    "Aromatic",
    "Hydrophobe",
    "LumpedHydrophobe",
)


@dataclass(frozen=True)
class SuCOSScore:
    score: float
    shape: float
    feature: float


_FDEF = AllChem.BuildFeatureFactory(
    os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
)


def _clamp_01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def _raw_features(mol: Chem.Mol):
    return [feat for feat in _FDEF.GetFeaturesForMol(mol) if feat.GetFamily() in KEEP_FEATURES]


def _feature_score(ref_mol: Chem.Mol, query_mol: Chem.Mol) -> float:
    ref_feats = _raw_features(ref_mol)
    query_feats = _raw_features(query_mol)

    if not ref_feats or not query_feats:
        return 0.0

    params = {name: FeatMaps.FeatMapParams() for name in KEEP_FEATURES}
    fmap = FeatMaps.FeatMap(
        feats=ref_feats,
        weights=[1.0] * len(ref_feats),
        params=params,
    )
    fmap.scoreMode = FeatMaps.FeatMapScoreMode.Best

    denom = float(min(len(ref_feats), len(query_feats)))
    if denom == 0.0:
        return 0.0

    return _clamp_01(fmap.ScoreFeats(query_feats) / denom)


def _shape_score(ref_mol: Chem.Mol, query_mol: Chem.Mol) -> float:
    protrude = rdShapeHelpers.ShapeProtrudeDist(
        ref_mol,
        query_mol,
        allowReordering=False,
    )
    return _clamp_01(1.0 - float(protrude))


def compute_sucos(ref_mol: Chem.Mol, query_mol: Chem.Mol) -> SuCOSScore:
    """Compute SuCOS score and components in [0,1]."""
    shape = _shape_score(ref_mol, query_mol)
    feature = _feature_score(ref_mol, query_mol)
    score = _clamp_01(0.5 * (shape + feature))
    return SuCOSScore(score=score, shape=shape, feature=feature)
