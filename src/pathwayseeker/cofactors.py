"""Cofactor and currency-metabolite set used by the oracle, reasoning and training layers.

This is the 42-compound set in the manuscript's Supplementary Information (cofactor
compound set table), used to generate the published training data
(training_data_generator_v3). The Methods text summarizes its main groups.
Graph construction in ``pathwayseeker.pipeline`` uses its own, smaller hub list.
"""

from typing import Dict, FrozenSet

COFACTOR_GROUPS: Dict[str, Dict[str, str]] = {
    "Water, ions and gases": {
        "C00001": "H2O", "C00007": "O2", "C00011": "CO2", "C00014": "NH3",
        "C00027": "H2O2", "C00080": "H+",
    },
    "Energy and phosphate carriers": {
        "C00002": "ATP", "C00008": "ADP", "C00020": "AMP", "C00075": "UTP",
        "C00009": "Orthophosphate", "C00013": "Diphosphate",
    },
    "Nucleotides": {
        "C00015": "UDP", "C00105": "UMP", "C00055": "CMP", "C00035": "GDP",
        "C00044": "GTP", "C00063": "CTP",
    },
    "Redox cofactors": {
        "C00003": "NAD+", "C00004": "NADH", "C00005": "NADPH", "C00006": "NADP+",
        "C00016": "FAD", "C01352": "FADH2", "C00061": "FMN",
    },
    "Acyl and one-carbon carriers": {
        "C00010": "CoA", "C00024": "Acetyl-CoA", "C00091": "Succinyl-CoA",
        "C00100": "Propanoyl-CoA", "C00101": "Tetrahydrofolate",
    },
    "Group donors": {
        "C00019": "S-Adenosyl-L-methionine", "C00021": "S-Adenosyl-L-homocysteine",
        "C00017": "Protein",
    },
    "High-degree metabolites": {
        "C00025": "L-Glutamate", "C00049": "L-Aspartate", "C00062": "L-Arginine",
        "C00065": "L-Serine",
    },
    "Sugars": {
        "C00029": "UDP-glucose", "C00031": "D-Glucose", "C00140": "N-Acetyl-D-glucosamine",
    },
    "Other": {
        "C00117": "D-Ribose 5-phosphate", "C00138": "Reduced ferredoxin",
    },
}

COFACTOR_NAMES: Dict[str, str] = {
    cid: name for group in COFACTOR_GROUPS.values() for cid, name in group.items()
}
COFACTORS: FrozenSet[str] = frozenset(COFACTOR_NAMES)

assert len(COFACTORS) == 42


def is_cofactor(compound_id: str) -> bool:
    return compound_id in COFACTORS
