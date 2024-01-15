from z3 import And, Or, sat, Solver
from comgen import SpeciesCollection
from comgen.constraint_system import TargetComposition, EMD, Synthesis, ONNX
from comgen.query import PETTIFOR_KEYS, element_to_pettifor, Query, get_radii
import pymatgen.core as pg
import csv

class IonicComposition(Query):
    def __init__(self, sps):
        super().__init__()
        assert isinstance(sps, SpeciesCollection)
        self._sps = sps
        self.precision = 0.1
        self._elmd_calculator = None
        self._setup()

    @property
    def precision(self):
        return self._precision

    @precision.setter
    def precision(self, val):
        self._precision = val

    @property
    def species(self):
        return self._sps

    def _setup(self):
        self.new_comp = TargetComposition(self._sps, self.constraints)
        self.new_comp.balance_charges()
        self.new_comp.restrict_charge_by_electronegativity()

    def _get_elmd_calc(self):
        if self._elmd_calculator is None:
            self._elmd_calculator = EMD(element_to_pettifor, PETTIFOR_KEYS, self.constraints)
        return self._elmd_calculator

    def elmd_close_to_one(self, compositions, bounds):
        if isinstance(bounds, float) or isinstance(bounds, int): [bounds] = bounds*len(compositions)
        assert len(bounds) == len(compositions)

        distances = []
        for comp, dist in zip(compositions, bounds):
            distances.append(self.new_comp.bound_distance(comp, self._get_elmd_calc(), ub=dist, return_constraint=True))
        self.constraints.append(Or(distances))    

    def elmd_far_from_all(self, compositions, bounds):
        if isinstance(bounds, float) or isinstance(bounds, int): bounds = [bounds]*len(compositions)
        assert len(bounds) == len(compositions)
        
        for comp, dist in zip(compositions, bounds):
            self.new_comp.bound_distance(comp, self._get_elmd_calc(), lb=dist)

    def made_from(self, ingredients):
        synthesis = Synthesis(ingredients, self.constraints)
        self.new_comp.synthesise_from(synthesis)

    def include_elements(self, elements, exact=None, *, lb=None, ub=None):
        self.new_comp.count_elements_from(elements, exact, lb=lb, ub=ub)

    def include_elements_quantity(self, elements, exact=None, *, lb=None, ub=None):
        self.new_comp.bound_elements_quantity(elements, exact, lb=lb, ub=ub)
    
    def distinct_elements(self, exact=None, *, lb=None, ub=None):
        self.new_comp.count_elements(exact, lb=lb, ub=ub)

    def total_atoms(self, exact=None, *, lb=None, ub=None):
        self.new_comp.count_atoms(exact, lb=lb, ub=ub)

    def exclude(self, compositions):
        for comp in compositions:
            self.new_comp.exclude_composition(comp)

    def radius_ratio(self, sps1, sps2, cn1=None, cn2=None, *, lb=None, ub=None):
        self.new_comp.include_species_pair_with_value_ratio(get_radii(sps1, cn1), get_radii(sps2, cn2), lb=lb, ub=ub)

    def category_prediction(self, onnx_model, category):
        model = ONNX(onnx_model, self.constraints)
        self.new_comp.property_predictor_category(model, category)

    def get_next(self, as_frac=False):
        model = super().get_next()
        if model is None:
            return None
        elt_quants = self.new_comp.format_solution(model, as_frac)
        self.new_comp.exclude_composition(elt_quants, self.precision)
        return {elt: str(q) for elt, q in elt_quants.items()}
