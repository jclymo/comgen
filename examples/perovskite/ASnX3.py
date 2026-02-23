from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import pymatgen.core as pg
from z3 import Int, Sum, And

from comgen import SpeciesCollection, PolyAtomicSpecies
from comgen.query.ionic import IonicComposition
from comgen.constraint_system.composition import UnitCell


@dataclass(frozen=True)
class Candidate:
    elements_frac: Dict[str, float]
    species_counts: Dict[str, int]
    score: float


def _float_from_model_val(v) -> float:
    # z3 rational -> float
    return float(v.numerator_as_long()) / float(v.denominator_as_long())


def build_species_space() -> SpeciesCollection:
    """
    SpeciesCollection can be built in multiple ways in comgen.
    This version is explicit and avoids relying on the packaged poly-ion file.
    """
    # Monatomic ions
    sn2 = pg.Species("Sn", 2)
    i_  = pg.Species("I", -1)
    br_ = pg.Species("Br", -1)
    cl_ = pg.Species("Cl", -1)

    # A-site options
    cs1 = pg.Species("Cs", 1)
    rb1 = pg.Species("Rb", 1)

    # Methylammonium: CH3NH3+ -> overall +1
    # Pymatgen composition uses element counts; charge stored in PolyAtomicSpecies( ..., oxi_state=+1 )
    ma_comp = pg.Composition({"C": 1, "H": 6, "N": 1})
    ma1 = PolyAtomicSpecies(ma_comp, 1)

    # Optional: formamidinium FA+ (CH5N2)+, common in perovskites
    fa_comp = pg.Composition({"C": 1, "H": 5, "N": 2})
    fa1 = PolyAtomicSpecies(fa_comp, 1)

    # Build the allowed species set
    return SpeciesCollection({sn2, i_, br_, cl_, cs1, rb1, ma1, fa1})


def add_abx3_constraints(query: IonicComposition, species: SpeciesCollection) -> None:
    """
    Enforce ABX3 in integer species counts:
      Sum(A) = k
      Sn2 = k
      Sum(X) = 3k
    and allow k = 1..4 formula units per cell.

    This is stricter and more chemically meaningful than only bounding fractional quantities.
    """
    # Create and attach a UnitCell to tie fractional quantities to integer counts.
    cell = UnitCell(species, query.constraints, query.return_vars)

    # Broad atom-count bounds so mixed organic/inorganic A-site still possible.
    # Needed because UnitCell.fit_composition requires bounds.  [oai_citation:4‡GitHub](https://raw.githubusercontent.com/jclymo/comgen/main/comgen/constraint_system/composition.py)
    cell.bound_total_atoms_count(lb=5, ub=80)
    query.target.fit_to_cell(cell)

    # Identify species names as strings (comgen stores vars keyed by str(sp))
    sn2 = str(pg.Species("Sn", 2))
    anions = {str(pg.Species("I", -1)), str(pg.Species("Br", -1)), str(pg.Species("Cl", -1))}
    a_cations = {str(pg.Species("Cs", 1)), str(pg.Species("Rb", 1)),
                 str(PolyAtomicSpecies(pg.Composition({"C": 1, "H": 6, "N": 1}), 1)),
                 str(PolyAtomicSpecies(pg.Composition({"C": 1, "H": 5, "N": 2}), 1))}

    k = Int("formula_units_k")
    query.constraints.append(And(k >= 1, k <= 4))

    # Pull the integer count vars
    sn_count = cell.species_count_vars(sn2)
    a_counts = [cell.species_count_vars(sp) for sp in a_cations]
    x_counts = [cell.species_count_vars(sp) for sp in anions]

    query.constraints.append(sn_count == k)
    query.constraints.append(Sum(a_counts) == k)
    query.constraints.append(Sum(x_counts) == 3 * k)

    # Optional: enforce that at least some A-site is organic (similarity to MA-based perovskites)
    ma_str = str(PolyAtomicSpecies(pg.Composition({"C": 1, "H": 6, "N": 1}), 1))
    query.constraints.append(cell.species_count_vars(ma_str) >= 0)  # keep allowed
    # If strict “must contain MA” is desired, uncomment:
    # query.constraints.append(cell.species_count_vars(ma_str) >= 1)


def stability_score(elements_frac: Dict[str, float]) -> float:
    """
    Proxy score: higher is "more stable".
    This is not thermodynamics. It is a triage heuristic.

    - Penalize large MA/FA organic content indirectly via C/H/N presence.
    - Reward some Br (often improves stability vs pure iodide).
    - Penalize Cl if large (often additive-level in many syntheses).

    Harsh reality: without an energy model or known structure, "stability" is a guess.
    """
    c = elements_frac.get("C", 0.0)
    h = elements_frac.get("H", 0.0)
    n = elements_frac.get("N", 0.0)

    i_ = elements_frac.get("I", 0.0)
    br = elements_frac.get("Br", 0.0)
    cl = elements_frac.get("Cl", 0.0)

    organic_penalty = 2.5 * (c + n) + 0.5 * h
    # Prefer halides dominated by I/Br with modest Br fraction
    # peak reward around br ~ 0.08-0.18 (tweakable)
    br_reward = 3.0 * (br - (br - 0.12) ** 2)
    cl_penalty = 4.0 * cl

    # Keep close to iodide baseline but allow some Br
    iodide_anchor = 1.5 * i_

    return iodide_anchor + br_reward - organic_penalty - cl_penalty


def generate_candidates(n: int = 100) -> List[Candidate]:
    species = build_species_space()
    q = IonicComposition(species)

    # Hard constraint: replace Pb with Sn, so disallow Pb entirely by construction.
    # q already includes charge balance and electronegativity heuristics internally.  [oai_citation:5‡GitHub](https://raw.githubusercontent.com/jclymo/comgen/main/comgen/constraint_system/composition.py)

    add_abx3_constraints(q, species)

    out: List[Candidate] = []
    for _ in range(n):
        model, _ = q.get_next()
        if model is None:
            break

        # element fractions (normalized)
        elt_frac = q.target.format_solution(model, as_frac=False)
        # pull integer species counts out of model for debugging / inspection
        species_counts = {}
        for v in q.return_vars:
            name = str(v)
            if "_speciescount" in name:
                mv = model[v]
                if mv is not None:
                    species_counts[name] = int(mv.as_long())

        out.append(Candidate(
            elements_frac=elt_frac,
            species_counts=species_counts,
            score=stability_score(elt_frac),
        ))

    out.sort(key=lambda c: c.score, reverse=True)
    return out


if __name__ == "__main__":
    cands = generate_candidates(n=200)
    for i, c in enumerate(cands[:20], start=1):
        print(f"\n#{i}  score={c.score:.3f}")
        print("elements:", c.elements_frac)