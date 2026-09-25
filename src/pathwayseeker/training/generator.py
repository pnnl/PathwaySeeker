#!/usr/bin/env python3
"""Generate fine-tuning examples from a graph.

This generator produced the released training set (data/training/training_v3.jsonl.gz),
ported from the unpublished research script training_data_generator_v3.py. It samples
examples of five evidence types (GRAPH_FACT, GRAPH_PATH, HYPOTHESIS, NO_PATH, INVALID) from
the enzyme-reaction-compound graph and applies the three-tier cofactor constraint defined by
L. M. O. Monteiro. With --balanced, NO_PATH and INVALID examples are capped at 20% of the
total.

Usage:
    pathwayseeker train-data --graph tversicolor --balanced --output train.jsonl

The prompt strings below match the released training data and should not be changed.
"""

from __future__ import annotations
import json
import random
import hashlib
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import networkx as nx

from pathwayseeker.cofactors import COFACTORS
from pathwayseeker.graph.multilayer import build_core_graph


class EvidenceType(str, Enum):
    """Evidence basis for every claim."""
    GRAPH_FACT = "GRAPH_FACT"
    GRAPH_PATH = "GRAPH_PATH"
    NO_PATH = "NO_PATH"
    HYPOTHESIS = "HYPOTHESIS"
    INVALID = "INVALID"


# =============================================================================
# COFACTOR DEFINITIONS
# =============================================================================

COMMON_COMPOUNDS: Set[str] = set(COFACTORS)  # 42-compound set, see pathwayseeker.cofactors

DEFAULT_COFACTOR_PENALTY = 10.0


def is_cofactor(compound_id: str) -> bool:
    """Check if a compound is a common cofactor/currency metabolite."""
    return compound_id in COMMON_COMPOUNDS


def get_cofactor_name(compound_id: str) -> str:
    """Get human-readable name for known cofactors."""
    COFACTOR_NAMES = {
        "C00001": "H2O", "C00007": "O2", "C00011": "CO2", "C00014": "NH3",
        "C00027": "H2O2", "C00080": "H+", "C00002": "ATP", "C00008": "ADP",
        "C00020": "AMP", "C00009": "Pi", "C00013": "PPi", "C00003": "NAD+",
        "C00004": "NADH", "C00005": "NADPH", "C00006": "NADP+", "C00016": "FAD",
        "C01352": "FADH2", "C00010": "CoA", "C00024": "Acetyl-CoA",
        "C00019": "SAM", "C00021": "SAH", "C00025": "Glutamate",
    }
    return COFACTOR_NAMES.get(compound_id, compound_id)


# =============================================================================
# SYSTEM PROMPT
# =============================================================================

SYSTEM_PROMPT = """You are a graph-constrained metabolic pathway assistant for multi-omics analysis.

GRAPH SCHEMA:
- COMPOUND: Metabolites identified by KEGG C-numbers (e.g., C00079 = L-Phenylalanine)
- REACTION: Biochemical transformations identified by KEGG R-numbers (e.g., R00686)
- ENZYME: Proteins identified by KEGG Orthology K-numbers (e.g., K00500)

EDGE TYPES:
- enzyme --[catalyzes]--> reaction
- reaction --[produces]--> compound  
- compound --[consumed_by]--> reaction

METABOLIC FLOW: A valid pathway requires substrate→product chains where each reaction's product serves as the next reaction's substrate.

COFACTOR POLICY:
Cofactors and currency metabolites (ATP, NAD+, CoA, H2O, etc.) require special handling:
1. ENDPOINT CONSTRAINT: Cofactors CANNOT be pathway start or end points. They do not define biological objectives.
2. INFERENCE CONSTRAINT: Inferred/hypothesized edges must NOT be cofactor-mediated. Only verified graph edges may involve cofactors.
3. TRAVERSAL PENALTY: Paths through cofactors are penalized. Prefer transformations of core metabolites (pathway backbone) over energy/redox currency exchanges.

Common cofactors include: H2O, O2, CO2, ATP/ADP/AMP, NAD+/NADH, NADP+/NADPH, FAD/FADH2, CoA and derivatives, SAM/SAH.

EVIDENCE TYPES (you must always specify one):
- GRAPH_FACT: Directly present in the provided graph context
- GRAPH_PATH: Connected via verified metabolic flow in the graph
- HYPOTHESIS: Proposed connection not present in graph (clearly marked as inference, cofactor-free)
- INVALID: The queried claim is false (including cofactor endpoint violations)

OUTPUT: Always respond with valid JSON containing evidence_type and confidence (0.0-1.0).
For HYPOTHESIS responses, explicitly separate verified_edges from inferred_edges, and ensure inferred_edges are not cofactor-mediated."""


# =============================================================================
# GRAPH INDEXING
# =============================================================================

@dataclass
class GraphIndexes:
    """Precomputed indexes for efficient sampling."""
    enzymes: Set[str] = field(default_factory=set)
    reactions: Set[str] = field(default_factory=set)
    compounds: Set[str] = field(default_factory=set)
    
    cofactors: Set[str] = field(default_factory=set)
    non_cofactor_compounds: Set[str] = field(default_factory=set)
    
    rxn_substrates: Dict[str, Set[str]] = field(default_factory=dict)
    rxn_products: Dict[str, Set[str]] = field(default_factory=dict)
    rxn_enzymes: Dict[str, Set[str]] = field(default_factory=dict)
    
    rxn_substrates_no_cofactors: Dict[str, Set[str]] = field(default_factory=dict)
    rxn_products_no_cofactors: Dict[str, Set[str]] = field(default_factory=dict)
    
    compound_consumed_by: Dict[str, Set[str]] = field(default_factory=dict)
    compound_produced_by: Dict[str, Set[str]] = field(default_factory=dict)
    
    enzyme_catalyzes: Dict[str, Set[str]] = field(default_factory=dict)
    
    components: List[Set[str]] = field(default_factory=list)
    node_to_component: Dict[str, int] = field(default_factory=dict)
    
    rxn_cofactor_participation: Dict[str, Set[str]] = field(default_factory=dict)


def build_indexes(G: nx.DiGraph, cofactor_aware: bool = True) -> GraphIndexes:
    """Build semantic indexes from the 3-layer graph."""
    idx = GraphIndexes()
    
    for node, data in G.nodes(data=True):
        node_type = data.get("type", "")
        if node_type == "enzyme":
            idx.enzymes.add(node)
            idx.enzyme_catalyzes[node] = set()
        elif node_type == "reaction":
            idx.reactions.add(node)
            idx.rxn_substrates[node] = set()
            idx.rxn_products[node] = set()
            idx.rxn_enzymes[node] = set()
            idx.rxn_substrates_no_cofactors[node] = set()
            idx.rxn_products_no_cofactors[node] = set()
            idx.rxn_cofactor_participation[node] = set()
        elif node_type == "compound":
            idx.compounds.add(node)
            idx.compound_consumed_by[node] = set()
            idx.compound_produced_by[node] = set()
            
            if cofactor_aware and is_cofactor(node):
                idx.cofactors.add(node)
            else:
                idx.non_cofactor_compounds.add(node)
    
    for u, v, data in G.edges(data=True):
        edge_type = data.get("type", "")
        
        if edge_type == "catalyzes" or edge_type == "catalyses":  # Handle both spellings
            if u in idx.enzymes and v in idx.reactions:
                idx.enzyme_catalyzes[u].add(v)
                idx.rxn_enzymes[v].add(u)
                
        elif edge_type == "produces":
            if u in idx.reactions and v in idx.compounds:
                idx.rxn_products[u].add(v)
                idx.compound_produced_by[v].add(u)
                
                if v not in idx.cofactors:
                    idx.rxn_products_no_cofactors[u].add(v)
                else:
                    idx.rxn_cofactor_participation[u].add(v)
                
        elif edge_type == "consumed_by":
            if u in idx.compounds and v in idx.reactions:
                idx.rxn_substrates[v].add(u)
                idx.compound_consumed_by[u].add(v)
                
                if u not in idx.cofactors:
                    idx.rxn_substrates_no_cofactors[v].add(u)
                else:
                    idx.rxn_cofactor_participation[v].add(u)
    
    undirected = G.to_undirected()
    idx.components = list(nx.connected_components(undirected))
    for i, component in enumerate(idx.components):
        for node in component:
            idx.node_to_component[node] = i
    
    return idx


# =============================================================================
# SERIALIZATION
# =============================================================================

def serialize_reaction(G: nx.DiGraph, rxn_id: str, idx: GraphIndexes, 
                       show_cofactor_participation: bool = True) -> str:
    """Serialize a reaction as a transformation unit."""
    substrates = idx.rxn_substrates.get(rxn_id, set())
    products = idx.rxn_products.get(rxn_id, set())
    enzymes = idx.rxn_enzymes.get(rxn_id, set())
    
    sub_strs = [f"{c} ({G.nodes[c].get('label', c)})" for c in sorted(substrates)]
    prod_strs = [f"{c} ({G.nodes[c].get('label', c)})" for c in sorted(products)]
    enz_strs = [f"{e} ({G.nodes[e].get('label', e)})" for e in sorted(enzymes)]
    
    lines = [
        f"REACTION: {rxn_id}",
        f"  SUBSTRATES: {', '.join(sub_strs) if sub_strs else 'unknown'}",
        f"  PRODUCTS: {', '.join(prod_strs) if prod_strs else 'unknown'}",
        f"  ENZYMES: {', '.join(enz_strs) if enz_strs else 'unknown'}",
    ]
    
    if show_cofactor_participation:
        cofactors_in_rxn = idx.rxn_cofactor_participation.get(rxn_id, set())
        if cofactors_in_rxn:
            cofactor_strs = [f"{c} ({get_cofactor_name(c)})" for c in sorted(cofactors_in_rxn)]
            lines.append(f"  COFACTORS: {', '.join(cofactor_strs)}")
    
    return "\n".join(lines)


def serialize_compound_context(G: nx.DiGraph, compound_id: str, idx: GraphIndexes,
                                indicate_cofactor: bool = True) -> str:
    """Serialize a compound with its metabolic context."""
    label = G.nodes[compound_id].get("label", compound_id)
    
    produced_by = idx.compound_produced_by.get(compound_id, set())
    consumed_by = idx.compound_consumed_by.get(compound_id, set())
    
    lines = [f"COMPOUND: {compound_id} ({label})"]
    
    if indicate_cofactor and is_cofactor(compound_id):
        lines.append(f"  TYPE: COFACTOR (currency metabolite)")
    
    lines.extend([
        f"  PRODUCED_BY: {', '.join(sorted(produced_by)) if produced_by else 'none'}",
        f"  CONSUMED_BY: {', '.join(sorted(consumed_by)) if consumed_by else 'none'}",
    ])
    return "\n".join(lines)


def serialize_metabolic_chain(
    G: nx.DiGraph, 
    chain: List[Tuple[str, str, str]],
    idx: GraphIndexes,
    annotate_cofactors: bool = True
) -> str:
    """Serialize a metabolic chain with full context."""
    lines = ["METABOLIC_CHAIN:"]
    for i, (substrate, reaction, product) in enumerate(chain):
        sub_label = G.nodes[substrate].get("label", substrate)
        prod_label = G.nodes[product].get("label", product)
        enzymes = idx.rxn_enzymes.get(reaction, set())
        enz_str = ", ".join(sorted(enzymes)) if enzymes else "unknown"
        
        step_line = (
            f"  STEP {i+1}: {substrate} ({sub_label}) "
            f"--[{reaction} via {enz_str}]--> "
            f"{product} ({prod_label})"
        )
        
        if annotate_cofactors:
            cofactors_here = []
            if is_cofactor(substrate):
                cofactors_here.append(f"{substrate}=cofactor")
            if is_cofactor(product):
                cofactors_here.append(f"{product}=cofactor")
            if cofactors_here:
                step_line += f" [{', '.join(cofactors_here)}]"
        
        lines.append(step_line)
    return "\n".join(lines)


def serialize_enzyme_context(G: nx.DiGraph, enzyme_id: str, idx: GraphIndexes) -> str:
    """Serialize an enzyme with all reactions it catalyzes."""
    label = G.nodes[enzyme_id].get("label", enzyme_id)
    reactions = idx.enzyme_catalyzes.get(enzyme_id, set())
    
    lines = [f"ENZYME: {enzyme_id} ({label})", "  CATALYZES:"]
    for rxn in sorted(reactions):
        subs = idx.rxn_substrates.get(rxn, set())
        prods = idx.rxn_products.get(rxn, set())
        lines.append(f"    {rxn}: {', '.join(sorted(subs))} -> {', '.join(sorted(prods))}")
    
    return "\n".join(lines)


# =============================================================================
# METABOLIC CHAIN FINDING
# =============================================================================

def find_metabolic_chains(
    G: nx.DiGraph, 
    idx: GraphIndexes,
    start_compound: str,
    max_length: int = 4,
    max_chains: int = 10,
    cofactor_aware: bool = True,
    cofactor_penalty: float = DEFAULT_COFACTOR_PENALTY
) -> List[Tuple[List[Tuple[str, str, str]], float]]:
    """Find valid metabolic chains starting from a compound."""
    chains_with_scores = []
    
    def dfs(current: str, current_chain: List[Tuple[str, str, str]], 
            visited: Set[str], current_score: float):
        if len(current_chain) >= max_length or len(chains_with_scores) >= max_chains:
            return
        
        consuming_rxns = idx.compound_consumed_by.get(current, set())
        
        for rxn in consuming_rxns:
            if rxn in visited:
                continue
            
            if cofactor_aware:
                products = idx.rxn_products_no_cofactors.get(rxn, set())
            else:
                products = idx.rxn_products.get(rxn, set())
            
            for product in products:
                if product in visited or product == current:
                    continue
                
                step_score = cofactor_penalty if is_cofactor(product) else 1.0
                new_score = current_score + step_score
                
                new_chain = current_chain + [(current, rxn, product)]
                chains_with_scores.append((new_chain, new_score))
                
                new_visited = visited | {rxn, product}
                dfs(product, new_chain, new_visited, new_score)
    
    if start_compound in idx.compounds:
        dfs(start_compound, [], {start_compound}, 0.0)
    
    chains_with_scores.sort(key=lambda x: (x[1], -len(x[0])))
    return chains_with_scores[:max_chains]


def find_all_metabolic_chains(
    G: nx.DiGraph, 
    idx: GraphIndexes,
    max_total_chains: int = 10000,
    min_chain_length: int = 2,
    cofactor_aware: bool = True,
    cofactor_penalty: float = DEFAULT_COFACTOR_PENALTY
) -> List[Tuple[List[Tuple[str, str, str]], float]]:
    """Find chains from all valid starting compounds."""
    all_chains = []
    
    if cofactor_aware:
        start_compounds = list(idx.non_cofactor_compounds)
    else:
        start_compounds = list(idx.compounds)
    
    random.shuffle(start_compounds)
    
    for start in start_compounds:
        if len(all_chains) >= max_total_chains:
            break
        
        chains = find_metabolic_chains(
            G, idx, start, max_length=5, max_chains=20,
            cofactor_aware=cofactor_aware, cofactor_penalty=cofactor_penalty
        )
        
        valid = [(c, s) for c, s in chains if len(c) >= min_chain_length]
        all_chains.extend(valid)
    
    return all_chains[:max_total_chains]


# =============================================================================
# SAMPLING FUNCTIONS
# =============================================================================

def sample_reaction_facts(
    G: nx.DiGraph, 
    idx: GraphIndexes, 
    n: int,
    cofactor_aware: bool = True
) -> List[Dict]:
    """Sample reaction-level facts (GRAPH_FACT evidence type).
    
    Generates multiple query types per reaction for variety:
    1. What are the products of reaction R?
    2. What are the substrates of reaction R?
    3. What reactions consume compound C?
    4. What reactions produce compound C?
    5. What enzymes catalyze reaction R?
    """
    examples = []
    seen_hashes = set()
    
    reactions = list(idx.reactions)
    random.shuffle(reactions)
    
    for rxn in reactions:
        if len(examples) >= n:
            break
        
        rxn_ctx = serialize_reaction(G, rxn, idx, show_cofactor_participation=cofactor_aware)
        products = list(idx.rxn_products.get(rxn, []))
        substrates = list(idx.rxn_substrates.get(rxn, []))
        enzymes = list(idx.rxn_enzymes.get(rxn, []))
        
        # Query Type 1: Products of reaction
        ex1 = {
            "input": f"GRAPH_CONTEXT:\n{rxn_ctx}\n\nQUERY: What are the products of {rxn}?",
            "output": json.dumps({
                "evidence_type": EvidenceType.GRAPH_FACT.value,
                "reaction": rxn,
                "products": products,
                "confidence": 1.0
            })
        }
        ex1_hash = hashlib.md5(ex1["input"].encode()).hexdigest()
        if ex1_hash not in seen_hashes:
            seen_hashes.add(ex1_hash)
            examples.append(ex1)
        
        # Query Type 2: Substrates of reaction
        if substrates and len(examples) < n:
            ex2 = {
                "input": f"GRAPH_CONTEXT:\n{rxn_ctx}\n\nQUERY: What are the substrates of {rxn}?",
                "output": json.dumps({
                    "evidence_type": EvidenceType.GRAPH_FACT.value,
                    "reaction": rxn,
                    "substrates": substrates,
                    "confidence": 1.0
                })
            }
            ex2_hash = hashlib.md5(ex2["input"].encode()).hexdigest()
            if ex2_hash not in seen_hashes:
                seen_hashes.add(ex2_hash)
                examples.append(ex2)
        
        # Query Type 3: Enzymes catalyzing reaction
        if enzymes and len(examples) < n:
            ex3 = {
                "input": f"GRAPH_CONTEXT:\n{rxn_ctx}\n\nQUERY: What enzymes catalyze {rxn}?",
                "output": json.dumps({
                    "evidence_type": EvidenceType.GRAPH_FACT.value,
                    "reaction": rxn,
                    "enzymes": enzymes,
                    "confidence": 1.0
                })
            }
            ex3_hash = hashlib.md5(ex3["input"].encode()).hexdigest()
            if ex3_hash not in seen_hashes:
                seen_hashes.add(ex3_hash)
                examples.append(ex3)
    
    # Query Types 4 & 5: Compound-centric queries (consumed_by, produced_by)
    compounds = list(idx.non_cofactor_compounds) if cofactor_aware else list(idx.compounds)
    random.shuffle(compounds)
    
    for compound in compounds:
        if len(examples) >= n:
            break
        
        comp_ctx = serialize_compound_context(G, compound, idx, indicate_cofactor=cofactor_aware)
        consumed_by = list(idx.compound_consumed_by.get(compound, []))
        produced_by = list(idx.compound_produced_by.get(compound, []))
        
        # Query Type 4: What reactions consume this compound?
        if consumed_by and len(examples) < n:
            ex4 = {
                "input": f"GRAPH_CONTEXT:\n{comp_ctx}\n\nQUERY: What reactions consume {compound}?",
                "output": json.dumps({
                    "evidence_type": EvidenceType.GRAPH_FACT.value,
                    "compound": compound,
                    "consumed_by_reactions": consumed_by,
                    "is_cofactor": is_cofactor(compound) if cofactor_aware else None,
                    "confidence": 1.0
                })
            }
            ex4_hash = hashlib.md5(ex4["input"].encode()).hexdigest()
            if ex4_hash not in seen_hashes:
                seen_hashes.add(ex4_hash)
                examples.append(ex4)
        
        # Query Type 5: What reactions produce this compound?
        if produced_by and len(examples) < n:
            ex5 = {
                "input": f"GRAPH_CONTEXT:\n{comp_ctx}\n\nQUERY: What reactions produce {compound}?",
                "output": json.dumps({
                    "evidence_type": EvidenceType.GRAPH_FACT.value,
                    "compound": compound,
                    "produced_by_reactions": produced_by,
                    "is_cofactor": is_cofactor(compound) if cofactor_aware else None,
                    "confidence": 1.0
                })
            }
            ex5_hash = hashlib.md5(ex5["input"].encode()).hexdigest()
            if ex5_hash not in seen_hashes:
                seen_hashes.add(ex5_hash)
                examples.append(ex5)
    
    return examples[:n]


def sample_metabolic_chains(
    G: nx.DiGraph, 
    idx: GraphIndexes, 
    n: int,
    cofactor_aware: bool = True,
    cofactor_penalty: float = DEFAULT_COFACTOR_PENALTY
) -> List[Dict]:
    """Sample valid metabolic chains (GRAPH_PATH evidence type)."""
    examples = []
    seen_hashes = set()
    
    compounds = list(idx.non_cofactor_compounds) if cofactor_aware else list(idx.compounds)
    random.shuffle(compounds)
    
    chains_collected = []
    for start in compounds:
        if len(chains_collected) >= n * 2:
            break
        
        chains = find_metabolic_chains(
            G, idx, start, max_length=5, max_chains=10,
            cofactor_aware=cofactor_aware, cofactor_penalty=cofactor_penalty
        )
        
        valid_chains = [(c, s) for c, s in chains if len(c) >= 2]
        chains_collected.extend(valid_chains)
    
    random.shuffle(chains_collected)
    
    for chain, score in chains_collected:
        if len(examples) >= n:
            break
        
        context = serialize_metabolic_chain(G, chain, idx, annotate_cofactors=cofactor_aware)
        
        start_compound = chain[0][0]
        end_compound = chain[-1][2]
        
        path = [start_compound]
        for sub, rxn, prod in chain:
            path.extend([rxn, prod])
        
        start_label = G.nodes[start_compound].get("label", start_compound)
        end_label = G.nodes[end_compound].get("label", end_compound)
        
        output_data = {
            "evidence_type": EvidenceType.GRAPH_PATH.value,
            "exists": True,
            "path": path,
            "num_steps": len(chain),
            "explanation": f"{start_compound} converts to {end_compound} via {len(chain)} enzymatic transformations.",
            "confidence": 1.0
        }
        
        if cofactor_aware:
            output_data["is_backbone_path"] = True
            output_data["path_score"] = score
        
        ex = {
            "input": f"GRAPH_CONTEXT:\n{context}\n\nQUERY: Is there a metabolic pathway from {start_compound} ({start_label}) to {end_compound} ({end_label})?",
            "output": json.dumps(output_data)
        }
        
        ex_hash = hashlib.md5(ex["input"].encode()).hexdigest()
        if ex_hash not in seen_hashes:
            seen_hashes.add(ex_hash)
            examples.append(ex)
    
    return examples[:n]


def sample_negative_pairs(
    G: nx.DiGraph, 
    idx: GraphIndexes, 
    n: int,
    cofactor_aware: bool = True
) -> List[Dict]:
    """Sample compound pairs with no metabolic connection (NO_PATH evidence type)."""
    examples = []
    seen_hashes = set()
    
    valid_compounds = list(idx.non_cofactor_compounds) if cofactor_aware else list(idx.compounds)
    
    compounds_by_component = {}
    for c in valid_compounds:
        comp_id = idx.node_to_component.get(c, -1)
        if comp_id not in compounds_by_component:
            compounds_by_component[comp_id] = []
        compounds_by_component[comp_id].append(c)
    
    component_ids = list(compounds_by_component.keys())
    
    if len(component_ids) < 2:
        return examples
    
    # Cross-component pairs (guaranteed no path)
    attempts = 0
    while len(examples) < n and attempts < n * 20:
        attempts += 1
        
        c1_id, c2_id = random.sample(component_ids, 2)
        compounds_c1 = compounds_by_component[c1_id]
        compounds_c2 = compounds_by_component[c2_id]
        
        if not compounds_c1 or not compounds_c2:
            continue
        
        src = random.choice(compounds_c1)
        tgt = random.choice(compounds_c2)
        
        pair_hash = hash(tuple(sorted([src, tgt])))
        if pair_hash in seen_hashes:
            continue
        seen_hashes.add(pair_hash)
        
        src_ctx = serialize_compound_context(G, src, idx, indicate_cofactor=cofactor_aware)
        tgt_ctx = serialize_compound_context(G, tgt, idx, indicate_cofactor=cofactor_aware)
        
        src_label = G.nodes[src].get("label", src)
        tgt_label = G.nodes[tgt].get("label", tgt)
        
        examples.append({
            "input": f"GRAPH_CONTEXT:\n{src_ctx}\n{tgt_ctx}\n\nQUERY: Is there a metabolic pathway from {src} ({src_label}) to {tgt} ({tgt_label})?",
            "output": json.dumps({
                "evidence_type": EvidenceType.NO_PATH.value,
                "exists": False,
                "reason": "Compounds are in disconnected components of the metabolic network.",
                "confidence": 1.0
            })
        })
    
    return examples[:n]


def sample_corrupted_chains(
    G: nx.DiGraph, 
    idx: GraphIndexes, 
    n: int,
    cofactor_aware: bool = True,
    cofactor_penalty: float = DEFAULT_COFACTOR_PENALTY
) -> List[Dict]:
    """Sample corrupted chains for hallucination detection (INVALID evidence type)."""
    examples = []
    seen_hashes = set()
    
    compounds = list(idx.non_cofactor_compounds) if cofactor_aware else list(idx.compounds)
    
    # First collect valid chains to corrupt
    valid_chains = []
    for _ in range(n * 3):
        start = random.choice(compounds)
        chains = find_metabolic_chains(
            G, idx, start, max_length=4, max_chains=5,
            cofactor_aware=cofactor_aware, cofactor_penalty=cofactor_penalty
        )
        valid_chains.extend([c for c, _ in chains if len(c) >= 2])
        if len(valid_chains) >= n * 2:
            break
    
    for chain in valid_chains:
        if len(examples) >= n:
            break
        
        corrupted = list(chain)
        error_idx = random.randint(0, len(chain) - 1)
        
        sub, rxn, prod = corrupted[error_idx]
        
        # Swap in wrong compound
        wrong_candidates = [c for c in compounds if c != sub and c != prod]
        if not wrong_candidates:
            continue
        
        wrong_compound = random.choice(wrong_candidates)
        if random.random() < 0.5:
            corrupted[error_idx] = (wrong_compound, rxn, prod)
        else:
            corrupted[error_idx] = (sub, rxn, wrong_compound)
        
        context = serialize_metabolic_chain(G, corrupted, idx, annotate_cofactors=cofactor_aware)
        
        start_compound = corrupted[0][0]
        end_compound = corrupted[-1][2]
        start_label = G.nodes.get(start_compound, {}).get("label", start_compound)
        end_label = G.nodes.get(end_compound, {}).get("label", end_compound)
        
        ex = {
            "input": f"GRAPH_CONTEXT:\n{context}\n\nQUERY: Verify this pathway from {start_compound} ({start_label}) to {end_compound} ({end_label}).",
            "output": json.dumps({
                "evidence_type": EvidenceType.INVALID.value,
                "valid": False,
                "error_at_step": error_idx + 1,
                "reason": "One or more steps in this pathway are not supported by the graph.",
                "confidence": 1.0
            })
        }
        
        ex_hash = hashlib.md5(ex["input"].encode()).hexdigest()
        if ex_hash not in seen_hashes:
            seen_hashes.add(ex_hash)
            examples.append(ex)
    
    return examples[:n]


def sample_enzyme_queries(
    G: nx.DiGraph, 
    idx: GraphIndexes, 
    n: int,
    cofactor_aware: bool = True
) -> List[Dict]:
    """Sample enzyme-centered queries."""
    examples = []
    seen_hashes = set()
    
    enzymes_with_reactions = [e for e in idx.enzymes if idx.enzyme_catalyzes.get(e)]
    random.shuffle(enzymes_with_reactions)
    
    for enzyme in enzymes_with_reactions:
        if len(examples) >= n:
            break
        
        ctx = serialize_enzyme_context(G, enzyme, idx)
        reactions = list(idx.enzyme_catalyzes.get(enzyme, []))
        
        ex = {
            "input": f"GRAPH_CONTEXT:\n{ctx}\n\nQUERY: What reactions does {enzyme} catalyze?",
            "output": json.dumps({
                "evidence_type": EvidenceType.GRAPH_FACT.value,
                "enzyme": enzyme,
                "catalyzes": reactions,
                "confidence": 1.0
            })
        }
        
        ex_hash = hashlib.md5(ex["input"].encode()).hexdigest()
        if ex_hash not in seen_hashes:
            seen_hashes.add(ex_hash)
            examples.append(ex)
    
    return examples[:n]


def sample_cofactor_endpoint_rejections(
    G: nx.DiGraph, 
    idx: GraphIndexes, 
    n: int
) -> List[Dict]:
    """Sample examples where cofactors are used as endpoints (INVALID)."""
    examples = []
    seen_hashes = set()
    
    cofactors = list(idx.cofactors)
    non_cofactors = list(idx.non_cofactor_compounds)
    
    if not cofactors or not non_cofactors:
        return examples
    
    for _ in range(n * 2):
        if len(examples) >= n:
            break
        
        cofactor = random.choice(cofactors)
        non_cofactor = random.choice(non_cofactors)
        
        # Cofactor as start
        if random.random() < 0.5:
            src, tgt = cofactor, non_cofactor
            reason = f"{src} ({get_cofactor_name(src)}) is a cofactor and cannot be a pathway starting point."
        else:
            src, tgt = non_cofactor, cofactor
            reason = f"{tgt} ({get_cofactor_name(tgt)}) is a cofactor and cannot be a pathway endpoint."
        
        src_ctx = serialize_compound_context(G, src, idx, indicate_cofactor=True)
        tgt_ctx = serialize_compound_context(G, tgt, idx, indicate_cofactor=True)
        
        src_label = G.nodes[src].get("label", src)
        tgt_label = G.nodes[tgt].get("label", tgt)
        
        ex = {
            "input": f"GRAPH_CONTEXT:\n{src_ctx}\n{tgt_ctx}\n\nQUERY: Is there a metabolic pathway from {src} ({src_label}) to {tgt} ({tgt_label})?",
            "output": json.dumps({
                "evidence_type": EvidenceType.INVALID.value,
                "exists": False,
                "reason": reason,
                "cofactor_policy_violation": True,
                "confidence": 1.0
            })
        }
        
        ex_hash = hashlib.md5(ex["input"].encode()).hexdigest()
        if ex_hash not in seen_hashes:
            seen_hashes.add(ex_hash)
            examples.append(ex)
    
    return examples[:n]


def sample_hypothesis_examples(
    G: nx.DiGraph, 
    idx: GraphIndexes, 
    n: int,
    cofactor_aware: bool = True
) -> List[Dict]:
    """Sample partial path hypotheses (HYPOTHESIS evidence type)."""
    examples = []
    seen_hashes = set()
    
    compounds = list(idx.non_cofactor_compounds) if cofactor_aware else list(idx.compounds)
    
    for _ in range(n * 3):
        if len(examples) >= n:
            break
        
        # Find a partial path and extend hypothetically
        start = random.choice(compounds)
        chains = find_metabolic_chains(
            G, idx, start, max_length=2, max_chains=5,
            cofactor_aware=cofactor_aware
        )
        
        if not chains:
            continue
        
        chain, _ = chains[0]
        if len(chain) < 1:
            continue
        
        # Get the endpoint and hypothesize extension to a random compound
        endpoint = chain[-1][2]
        
        # Find a compound NOT reachable from endpoint
        target_candidates = [c for c in compounds if c != endpoint]
        if not target_candidates:
            continue
        
        target = random.choice(target_candidates)
        
        # Check if actually unreachable
        target_chains = find_metabolic_chains(G, idx, endpoint, max_length=4, max_chains=5)
        reachable = any(c[-1][2] == target for c, _ in target_chains)
        
        if reachable:
            continue  # Don't hypothesize about reachable compounds
        
        context = serialize_metabolic_chain(G, chain, idx, annotate_cofactors=cofactor_aware)
        
        start_label = G.nodes[start].get("label", start)
        target_label = G.nodes[target].get("label", target)
        
        verified_edges = []
        for sub, rxn, prod in chain:
            verified_edges.append({
                "from": sub, "to": prod, "reaction": rxn, "confidence": 1.0
            })
        
        ex = {
            "input": f"GRAPH_CONTEXT:\n{context}\n\nQUERY: Is there a metabolic pathway from {start} ({start_label}) to {target} ({target_label})?",
            "output": json.dumps({
                "evidence_type": EvidenceType.HYPOTHESIS.value,
                "exists": True,
                "verified_edges": verified_edges,
                "inferred_edges": [{
                    "from": endpoint,
                    "to": target,
                    "confidence": 0.6,
                    "basis": "Biochemically plausible transformation, but not verified in current graph."
                }],
                "overall_confidence": 0.6,
                "explanation": f"Partial path verified from {start} to {endpoint}. Extension to {target} is hypothesized."
            })
        }
        
        ex_hash = hashlib.md5(ex["input"].encode()).hexdigest()
        if ex_hash not in seen_hashes:
            seen_hashes.add(ex_hash)
            examples.append(ex)
    
    return examples[:n]


def sample_metabolic_chains_exhaustive(
    G: nx.DiGraph, 
    idx: GraphIndexes, 
    cofactor_aware: bool = True,
    cofactor_penalty: float = DEFAULT_COFACTOR_PENALTY,
    max_total: int = 15000
) -> List[Dict]:
    """
    Exhaustively extract metabolic chains from ALL starting compounds.
    This maximizes GRAPH_PATH examples from the graph.
    """
    examples = []
    seen_hashes = set()
    seen_pairs = set()  # Track (src, tgt) to avoid duplicate pathways
    
    compounds = list(idx.non_cofactor_compounds) if cofactor_aware else list(idx.compounds)
    random.shuffle(compounds)
    
    print(f"    Searching from {len(compounds)} compounds...")
    
    chains_found = 0
    for i, start in enumerate(compounds):
        if len(examples) >= max_total:
            break
        
        # Find chains from this starting compound
        chains = find_metabolic_chains(
            G, idx, start, 
            max_length=6,       # Longer chains
            max_chains=50,      # More chains per start
            cofactor_aware=cofactor_aware, 
            cofactor_penalty=cofactor_penalty
        )
        
        for chain, score in chains:
            if len(chain) < 2:
                continue
            
            start_compound = chain[0][0]
            end_compound = chain[-1][2]
            
            # Skip if we've seen this pair
            pair = (start_compound, end_compound)
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            
            context = serialize_metabolic_chain(G, chain, idx, annotate_cofactors=cofactor_aware)
            
            path = [start_compound]
            for sub, rxn, prod in chain:
                path.extend([rxn, prod])
            
            start_label = G.nodes[start_compound].get("label", start_compound)
            end_label = G.nodes[end_compound].get("label", end_compound)
            
            output_data = {
                "evidence_type": EvidenceType.GRAPH_PATH.value,
                "exists": True,
                "path": path,
                "num_steps": len(chain),
                "explanation": f"{start_compound} converts to {end_compound} via {len(chain)} enzymatic transformations.",
                "confidence": 1.0
            }
            
            if cofactor_aware:
                output_data["is_backbone_path"] = True
                output_data["path_score"] = score
            
            ex = {
                "input": f"GRAPH_CONTEXT:\n{context}\n\nQUERY: Is there a metabolic pathway from {start_compound} ({start_label}) to {end_compound} ({end_label})?",
                "output": json.dumps(output_data)
            }
            
            ex_hash = hashlib.md5(ex["input"].encode()).hexdigest()
            if ex_hash not in seen_hashes:
                seen_hashes.add(ex_hash)
                examples.append(ex)
                chains_found += 1
        
        # Progress update
        if (i + 1) % 200 == 0:
            print(f"    ... processed {i+1}/{len(compounds)} compounds, found {chains_found} chains")
    
    print(f"    Total chains extracted: {len(examples)}")
    return examples


# =============================================================================
# MAIN GENERATION WITH BALANCE ENFORCEMENT
# =============================================================================

def generate_training_data(
    G: nx.DiGraph,
    n_total: int = 30000,
    max_negative_ratio: float = 0.25,
    cofactor_aware: bool = True,
    cofactor_penalty: float = DEFAULT_COFACTOR_PENALTY,
    seed: int = 42
) -> Tuple[List[Dict], Dict]:
    """
    Generate balanced training dataset with post-hoc rebalancing.
    
    Strategy: Extract ALL possible positive examples first, then cap negatives.
    This ensures we maximize the graph's potential for positive training signal.
    
    Args:
        G: Multi-omics graph
        n_total: Target total examples (will extract up to this, then rebalance)
        max_negative_ratio: Maximum ratio for NO_PATH + INVALID combined (default 0.25)
        cofactor_aware: Enable cofactor constraints
        cofactor_penalty: Penalty for cofactor traversal
        seed: Random seed
    
    Returns:
        (examples, stats)
    """
    random.seed(seed)
    
    idx = build_indexes(G, cofactor_aware=cofactor_aware)
    
    stats = {
        "graph": {
            "enzymes": len(idx.enzymes),
            "reactions": len(idx.reactions),
            "compounds": len(idx.compounds),
            "cofactors": len(idx.cofactors),
            "non_cofactor_compounds": len(idx.non_cofactor_compounds),
            "components": len(idx.components),
        },
        "config": {
            "n_total": n_total,
            "max_negative_ratio": max_negative_ratio,
            "cofactor_aware": cofactor_aware,
            "cofactor_penalty": cofactor_penalty,
            "seed": seed,
        },
        "samples": {}
    }
    
    print(f"Graph: {stats['graph']['enzymes']} enzymes, "
          f"{stats['graph']['reactions']} reactions, "
          f"{stats['graph']['compounds']} compounds "
          f"({stats['graph']['cofactors']} cofactors, {stats['graph']['non_cofactor_compounds']} backbone), "
          f"{stats['graph']['components']} components")
    print(f"\nStrategy: Extract maximum positives, then cap negatives at {max_negative_ratio:.0%}")
    
    # =========================================================================
    # PHASE 1: EXTRACT ALL POSITIVE EXAMPLES (maximize)
    # =========================================================================
    print(f"\n{'='*60}")
    print("PHASE 1: EXTRACTING POSITIVE EXAMPLES (MAXIMIZE)")
    print('='*60)
    
    all_examples = []
    
    # 1. Reaction Facts - extract ALL (this is our core knowledge)
    print(f"\n[1] Extracting ALL reaction facts...")
    reaction_facts = sample_reaction_facts(G, idx, n=20000, cofactor_aware=cofactor_aware)
    all_examples.extend(reaction_facts)
    stats["samples"]["reaction_facts"] = len(reaction_facts)
    print(f"  → Extracted {len(reaction_facts)} GRAPH_FACT examples")
    
    # 2. Metabolic Chains - extract aggressively
    print(f"\n[2] Extracting metabolic chains (aggressive)...")
    chains = sample_metabolic_chains_exhaustive(G, idx, cofactor_aware=cofactor_aware, cofactor_penalty=cofactor_penalty)
    all_examples.extend(chains)
    stats["samples"]["metabolic_chains"] = len(chains)
    print(f"  → Extracted {len(chains)} GRAPH_PATH examples")
    
    # 3. Enzyme Queries
    print(f"\n[3] Extracting enzyme queries...")
    enzymes_with_reactions = [e for e in idx.enzymes if idx.enzyme_catalyzes.get(e)]
    print(f"    Enzymes with catalyzes edges: {len(enzymes_with_reactions)} of {len(idx.enzymes)}")
    if len(enzymes_with_reactions) == 0:
        print(f"    Warning: no enzyme-catalyzes edges found. Check graph edge types.")
        # Debug: show some edges
        edge_types = {}
        for u, v, data in G.edges(data=True):
            et = data.get("type", "unknown")
            edge_types[et] = edge_types.get(et, 0) + 1
        print(f"    Edge types in graph: {edge_types}")
    enzyme_qs = sample_enzyme_queries(G, idx, n=5000, cofactor_aware=cofactor_aware)
    all_examples.extend(enzyme_qs)
    stats["samples"]["enzyme_queries"] = len(enzyme_qs)
    print(f"  → Extracted {len(enzyme_qs)} enzyme query examples")
    
    # 4. Hypothesis examples (CAPPED at 5% of expected total)
    # Don't let HYPOTHESIS dominate just because other categories are small
    max_hypothesis = int(len(all_examples) * 0.08)  # Cap at ~8% of what we have so far
    print(f"\n[4] Extracting hypothesis examples (capped at {max_hypothesis})...")
    hypotheses = sample_hypothesis_examples(G, idx, n=max_hypothesis, cofactor_aware=cofactor_aware)
    all_examples.extend(hypotheses)
    stats["samples"]["hypotheses"] = len(hypotheses)
    print(f"  → Extracted {len(hypotheses)} HYPOTHESIS examples")
    
    # Count positives so far
    positive_count = len(reaction_facts) + len(chains) + len(enzyme_qs) + len(hypotheses)
    print(f"\n  Total positive examples: {positive_count}")
    
    # =========================================================================
    # PHASE 2: ADD CAPPED NEGATIVE EXAMPLES
    # =========================================================================
    print(f"\n{'='*60}")
    print("PHASE 2: ADDING CAPPED NEGATIVE EXAMPLES")
    print('='*60)
    
    # Calculate max negatives based on current positives
    # If we want negatives to be max_negative_ratio of total:
    # negatives / (positives + negatives) = max_negative_ratio
    # negatives = positives * max_negative_ratio / (1 - max_negative_ratio)
    max_negatives = int(positive_count * max_negative_ratio / (1 - max_negative_ratio))
    
    # Split between NO_PATH and INVALID
    n_nopath = int(max_negatives * 0.55)  # 55% of negatives = NO_PATH
    n_invalid = max_negatives - n_nopath   # 45% of negatives = INVALID
    n_cofactor_rejections = int(n_invalid * 0.4) if cofactor_aware else 0
    n_corrupted = n_invalid - n_cofactor_rejections
    
    print(f"\n  Max negatives allowed: {max_negatives} ({max_negative_ratio:.0%} of total)")
    print(f"    - NO_PATH target: {n_nopath}")
    print(f"    - INVALID target: {n_invalid} (corrupted: {n_corrupted}, cofactor: {n_cofactor_rejections})")
    
    # 5. Negative Pairs (NO_PATH) - CAPPED
    print(f"\n[5] Sampling negative pairs (capped at {n_nopath})...")
    negatives = sample_negative_pairs(G, idx, n_nopath, cofactor_aware=cofactor_aware)
    all_examples.extend(negatives)
    stats["samples"]["negative_pairs"] = len(negatives)
    print(f"  → Generated {len(negatives)} NO_PATH examples")
    
    # 6. Corrupted Chains (INVALID) - CAPPED
    print(f"\n[6] Sampling corrupted chains (capped at {n_corrupted})...")
    corrupted = sample_corrupted_chains(G, idx, n_corrupted, cofactor_aware=cofactor_aware, cofactor_penalty=cofactor_penalty)
    all_examples.extend(corrupted)
    stats["samples"]["corrupted_chains"] = len(corrupted)
    print(f"  → Generated {len(corrupted)} INVALID (corrupted) examples")
    
    # 7. Cofactor Rejections (INVALID)
    if cofactor_aware and n_cofactor_rejections > 0:
        print(f"\n[7] Sampling cofactor rejections (capped at {n_cofactor_rejections})...")
        cofactor_rejections = sample_cofactor_endpoint_rejections(G, idx, n_cofactor_rejections)
        all_examples.extend(cofactor_rejections)
        stats["samples"]["cofactor_rejections"] = len(cofactor_rejections)
        print(f"  → Generated {len(cofactor_rejections)} INVALID (cofactor) examples")
    else:
        stats["samples"]["cofactor_rejections"] = 0
    
    # =========================================================================
    # PHASE 3: DEDUPLICATE AND CALCULATE FINAL DISTRIBUTION
    # =========================================================================
    print(f"\n{'='*60}")
    print("PHASE 3: DEDUPLICATION AND FINAL STATS")
    print('='*60)
    
    # Deduplicate
    seen = set()
    deduped = []
    for ex in all_examples:
        ex_hash = hashlib.md5(ex["input"].encode()).hexdigest()
        if ex_hash not in seen:
            seen.add(ex_hash)
            deduped.append(ex)
    
    stats["samples"]["before_dedup"] = len(all_examples)
    stats["samples"]["after_dedup"] = len(deduped)
    stats["samples"]["duplicates_removed"] = len(all_examples) - len(deduped)
    
    random.shuffle(deduped)
    
    # Calculate actual distribution
    evidence_counts = {}
    for ex in deduped:
        try:
            output = json.loads(ex["output"])
            ev_type = output.get("evidence_type", "UNKNOWN")
            evidence_counts[ev_type] = evidence_counts.get(ev_type, 0) + 1
        except:
            pass
    
    stats["samples"]["total"] = len(deduped)
    stats["samples"]["by_evidence_type"] = evidence_counts
    
    # Calculate positive/negative ratios
    total = len(deduped)
    positive = evidence_counts.get("GRAPH_FACT", 0) + evidence_counts.get("GRAPH_PATH", 0)
    negative = evidence_counts.get("NO_PATH", 0) + evidence_counts.get("INVALID", 0)
    
    stats["samples"]["positive_ratio"] = positive / total if total > 0 else 0
    stats["samples"]["negative_ratio"] = negative / total if total > 0 else 0
    
    print(f"\n{'='*60}")
    print("GENERATION SUMMARY")
    print('='*60)
    print(f"Total examples: {len(deduped)}")
    print(f"\nDistribution:")
    for ev_type, count in sorted(evidence_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / total if total > 0 else 0
        print(f"  {ev_type}: {count} ({pct:.1f}%)")
    print(f"\nPositive (GRAPH_FACT + GRAPH_PATH): {positive} ({100*positive/total:.1f}%)")
    print(f"Negative (NO_PATH + INVALID): {negative} ({100*negative/total:.1f}%)")
    
    # Warn if still imbalanced
    if negative / total > max_negative_ratio + 0.05:
        print(f"\nWarning: negative ratio ({100*negative/total:.1f}%) exceeds target ({100*max_negative_ratio:.0f}%)")
        print(f"    Graph may not support enough positive examples.")
    
    return deduped, stats


def convert_to_openai_format(
    examples: List[Dict], 
    system_prompt: str = SYSTEM_PROMPT
) -> List[Dict]:
    """Convert to OpenAI fine-tuning JSONL format."""
    return [{
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": ex["input"]},
            {"role": "assistant", "content": ex["output"]}
        ]
    } for ex in examples]


def save_training_data(
    examples: List[Dict],
    output_path: Path,
    system_prompt: str = SYSTEM_PROMPT
):
    """Save training data to JSONL file in OpenAI format."""
    openai_examples = convert_to_openai_format(examples, system_prompt)
    
    with open(output_path, "w") as f:
        for ex in openai_examples:
            f.write(json.dumps(ex) + "\n")
    
    print(f"Saved {len(openai_examples)} examples to {output_path}")


# =============================================================================
# CLI
# =============================================================================

def main(argv=None):
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate balanced training data for graph-constrained LLM fine-tuning (v3)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  pathwayseeker train-data --graph tversicolor --balanced --output train.jsonl
        """
    )
    
    parser.add_argument("--output", type=Path, default=Path("training_data_v3.jsonl"),
                        help="Output JSONL file")
    parser.add_argument("--stats-output", type=Path, default=None,
                        help="Output stats JSON file (default: <output>_stats.json)")
    
    # Balance controls
    parser.add_argument("--balanced", action="store_true",
                        help="Use balanced preset (recommended for production)")
    parser.add_argument("--n-total", type=int, default=30000,
                        help="Target total examples (positives extracted first, negatives capped)")
    parser.add_argument("--max-negative-ratio", type=float, default=0.25,
                        help="Maximum negative ratio (NO_PATH + INVALID), default 0.25")
    
    # Cofactor controls
    parser.add_argument("--no-cofactor-aware", action="store_true",
                        help="Disable cofactor awareness (for ablation)")
    parser.add_argument("--cofactor-penalty", type=float, default=DEFAULT_COFACTOR_PENALTY,
                        help=f"Penalty for traversing cofactor nodes (default: {DEFAULT_COFACTOR_PENALTY})")
    
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility")
    
    parser.add_argument("--graph-dir", type=Path, default=Path("data/output"),
                        help="Directory with the pipeline outputs that define the graph")

    args = parser.parse_args(argv)
    
    # Set stats output path
    if args.stats_output is None:
        args.stats_output = args.output.with_suffix('.stats.json')
    
    print("Building graph...")
    G = build_core_graph(str(args.graph_dir))
    
    cofactor_aware = not args.no_cofactor_aware
    
    # Balanced preset
    if args.balanced:
        args.max_negative_ratio = 0.20
    
    print(f"\nConfiguration:")
    print(f"  Max negative ratio: {args.max_negative_ratio:.0%}")
    print(f"  Cofactor-aware: {cofactor_aware}")
    print(f"  Cofactor penalty: {args.cofactor_penalty}")
    print(f"  Seed: {args.seed}")
    
    examples, stats = generate_training_data(
        G,
        n_total=args.n_total,
        max_negative_ratio=args.max_negative_ratio,
        cofactor_aware=cofactor_aware,
        cofactor_penalty=args.cofactor_penalty,
        seed=args.seed
    )
    
    save_training_data(examples, args.output, SYSTEM_PROMPT)
    
    with open(args.stats_output, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Stats saved to {args.stats_output}")


if __name__ == "__main__":
    main()
