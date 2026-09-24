# Analytics Module

Analysis and visualization tools.

::: cofolder.modules.analytics
    options:
      show_root_heading: true
      show_source: true
      heading_level: 2

## Typed interaction fingerprints

::: cofolder.modules.analytics.reference_ifp
    options:
      show_root_heading: true
      show_source: true
      heading_level: 3
      members:
        - AtomIdentity
        - InteractionGeometry
        - InteractionEvent
        - InteractionFingerprint
        - InteractionKey
        - IFPExtractionConfig
        - extract_interaction_fingerprint

`InteractionFingerprint.interactions` is the deduplicated residue/type feature set.
`InteractionFingerprint.events` contains every atom-level ProLIF occurrence. Event
geometry uses Å for distances and degrees for angles. The extractor accepts PDB and
mmCIF paths; mmCIF conversion remains contained in the isolated ProLIF worker.

```python
from pathlib import Path

from cofolder.modules.analytics.reference_ifp import (
    IFPExtractionConfig,
    IFPTaxonomy,
    LigandSelector,
    extract_interaction_fingerprint,
)

fingerprint = extract_interaction_fingerprint(
    Path("prediction.cif"),
    ligand=LigandSelector(chain_id="B"),
    receptor_chains=("A",),
    config=IFPExtractionConfig(taxonomy=IFPTaxonomy.PROLIF),
)
asp168_backbone_n = [
    event
    for event in fingerprint.events
    if any(
        atom.residue_name == "ASP"
        and atom.residue_number == 168
        and atom.atom_name == "N"
        for atom in event.protein_atoms
    )
]
```

## Interaction-fingerprint clustering

::: cofolder.modules.analytics.ifp_clustering
    options:
      show_root_heading: true
      show_source: true
      heading_level: 3
      members:
        - IFPClusteringResult
        - cluster_binary_ifps
        - cluster_interaction_fingerprints
