from z3 import And, Or, Not, Sum, Real, Q
from comgen import SpeciesCollection, PolyAtomicSpecies
from comgen.constraint_system.common import zero_weighted_sum, apply_bounds, check_bounds, bound_weighted_average_value_ratio
import pymatgen.core as pg
from fractions import Fraction

class TargetComposition:
    def __init__(self, permitted_species, constraint_log, return_vars):
        if not isinstance(permitted_species, SpeciesCollection): 
            raise TypeError("permitted_species argument must be a SpeciesCollection.")
        
        self.name = f'TargetComposition{str(id(self))}'
        
        self.permitted_species = permitted_species

        self.cons = constraint_log
        self.return_vars = return_vars

        self.element_quantity_variable_collection = {} 
        self.species_quantity_variable_collection = {}
        self._setup()

    def _new_element_quantity_var(self, el):
        el_id = str(el)
        var = Real(f'{self.name}_{el_id}_elementquantity')
        self.element_quantity_variable_collection[el_id] = var
        self.return_vars.append(var)

    def _new_species_quantity_var(self, sp):
        sp_id = str(sp)
        var = Real(f'{self.name}_{sp_id}_speciesquantity')
        self.species_quantity_variable_collection[sp_id] = var
        self.return_vars.append(var)

    def _setup(self):
        for sp in self.permitted_species:
            self._new_species_quantity_var(sp)
        for el in self.elements:
            self._new_element_quantity_var(el)
        
        # add basic constraints to give correct semantics to the objects variables
        # 1. element and species quantities must be non-negative
        for sp in self.permitted_species:
            self.cons.append(self.species_quantity_vars(sp) >= 0)
        for el in self.elements:
            self.cons.append(self.element_quantity_vars(el) >= 0)

        # 2. total species quantity equals corresponding element quantity
        for el, sps in self.permitted_species.group_by_element_view().items():
            sps_vars, sps_weights = [], []
            for sp in sps:
                sps_vars.append(self.species_quantity_vars(sp))
                weight = 1
                if isinstance(sp, PolyAtomicSpecies):
                    weight = sp.multiplier(el)
                sps_weights.append(weight)
            
            self.cons.append(Sum([var*weight for var, weight in zip(sps_vars, sps_weights)]) == self.element_quantity_vars(el))

        # 3. total element quantity is 1
        vars = [self.element_quantity_vars(el) for el in self.elements]
        self.cons.append(Sum(*vars) == 1)

    @property
    def elements(self):
        return self.permitted_species.group_by_element_view().keys()

    def element_quantity_vars(self, el=None): # ideally these wouldn't be available as public methods. with a reason we can pass out the vars but what they represent shouldn't be known. 
        if el is not None and not isinstance(el, str): el = str(el)

        if el:
            return self.element_quantity_variable_collection[el]
        return self.element_quantity_variable_collection
    
    def species_quantity_vars(self, sp=None):
        if sp is not None and not isinstance(sp, str): sp = str(sp)

        if sp:
            return self.species_quantity_variable_collection[sp]
        return self.species_quantity_variable_collection    

    def balance_charges(self, return_constraint=False):
        """Weighted sum of atom charges is zero.
            respect_electronegativity: if this is True then we require that 

        """
        sps_quants = [self.species_quantity_vars(str(sp)) for sp in self.permitted_species]
        sps_charges = [sp.oxi_state for sp in self.permitted_species]
        balance_charge_cons = zero_weighted_sum(sps_quants, sps_charges)

        if return_constraint:
            return balance_charge_cons 
        self.cons.append(balance_charge_cons)

    def restrict_charge_by_electronegativity(self, return_constraint=False):
        """ no element can have both positively and negatively charged ions in the same composition
            all positively charged ions must have lower e-neg than all negatively charged ions
        """
        elt_grouped_sps = self.permitted_species.group_by_element_view()
        eneg_cons = []
        
        pos, neg = {}, {}
        for elt, sps in elt_grouped_sps.items():
            sps = {sp for sp in sps if isinstance(sp, pg.Species)} # remove polyatomic species
            pos[str(elt)] = Or([self.species_quantity_vars(str(sp)) > 0 for sp in sps if sp.oxi_state > 0])
            neg[str(elt)] = Or([self.species_quantity_vars(str(sp)) > 0 for sp in sps if sp.oxi_state < 0])

        for el_1 in elt_grouped_sps.keys():
            for el_2 in elt_grouped_sps.keys():
                if el_1.X > el_2.X or el_1 == el_2:
                    eneg_cons.append(Not(And(pos[str(el_1)], neg[str(el_2)])))
                
        if return_constraint:
            return And(eneg_cons)
        self.cons.append(And(eneg_cons))

    def count_elements_from(self, elements: set, exact: int=None, return_constraint=False, *, lb: int=None, ub: int=None):
        """Constrain the number of elements included from the given sets 
        """
        count_elts_present = Sum([self.element_quantity_vars(str(elt)) > 0 for elt in elements])
        bound_elts_present_cons = apply_bounds(count_elts_present, exact, lb=lb, ub=ub)
        
        if return_constraint:
            return bound_elts_present_cons
        self.cons.append(bound_elts_present_cons)
    
    def count_elements(self, exact: int=None, return_constraint=False, *, lb: int=None, ub: int=None):
        """Constrain the number of elements included in the composition
        """
        elements = {str(elt) for elt in self.elements}
        return self.count_elements_from(elements, exact, return_constraint, lb=lb, ub=ub)

    def count_atoms(self, exact: int=None, return_constraint=False, *, lb: int=None, ub: int=None):
        """Constrain the total number of atoms. 
        The quantity of each species must be rational with denominator matching number of atoms. 
        """
        check_bounds(exact, lb, ub)
        
        sps_vars = self.species_quantity_vars()
        if exact: lb, ub = exact, exact
        if lb is None: lb = 1
        if ub is None: raise ValueError('Please provide an upper bound on the number of atoms.')

        select_for_total = []
        for num_atoms in range(lb, ub+1):
            select_n = []
            for var in sps_vars.values():
                select_n.append(Or([var == Q(n, num_atoms) for n in range(num_atoms+1)]))
            select_for_total.append(And(select_n))

        select_num_atoms_cons = Or(select_for_total)
        if return_constraint:
            return select_num_atoms_cons
        self.cons.append(select_num_atoms_cons)

    def bound_elements_quantity(self, elements: set, exact: float=None, return_constraint=False, *, lb: float=None, ub: float=None):
        """Constraint the total quantity across elements in the given set.
        """
        total_elts_quantity = Sum([self.element_quantity_vars(str(elt)) for elt in elements])
        bound_elts_quantity_cons = apply_bounds(total_elts_quantity, exact, lb=lb, ub=ub)

        if return_constraint:
            return bound_elts_quantity_cons
        self.cons.append(bound_elts_quantity_cons)

    def bound_species_quantity(self, sps: set, exact: float=None, return_constraint=False, *, lb: float=None, ub: float=None):
        """Constraint the total quantity across species in the given set.
        """
        total_sps_quantity = Sum([self.species_quantity_vars(str(sp)) for sp in sps])
        bound_sps_quantity_cons = apply_bounds(total_sps_quantity, exact, lb=lb, ub=ub)

        if return_constraint:
            return bound_sps_quantity_cons
        self.cons.append(bound_sps_quantity_cons)

    def exclude_composition(self, composition, precision=0.1, return_constraint=False):
        """Exclude a composition from the composition space. 
        
        Args:
            composition (dict): composition to exclude
            precision (float): tolerance for excluding composition
        """
        if isinstance(composition, pg.Composition):
            composition = dict(composition)
        
        cons = []
        for elt, quant in composition.items():
            lb, ub = quant - precision, quant + precision
            # get constraint that fixes quantity close to this solution
            cons.append(apply_bounds(self.element_quantity_vars(elt), lb=lb, ub=ub))
        
        exclude_cons = Not(And(cons)) # not all quantities are close to this solution
        if return_constraint:
            return exclude_cons
        self.cons.append(exclude_cons)

    def select_species_pair(self, pairs, return_constraint=False):
        sps_quants = self.species_quantity_vars()
        select_cons = []
        for sp1, sp2 in pairs:
            select_cons.append(And(sps_quants[sp1] > 0, sps_quants[sp2] > 0))
        
        if return_constraint:
            return Or(select_cons)
        self.cons.append(Or(select_cons))

    def bound_average_species_value_ratio(self, sps_1, sps_2, return_constraint=False, *, lb: float=None, ub: float=None):
        """Constrain the ratio of the average value of species in two sets.
        """
        sps_quants = self.species_quantity_vars()
        vars_1 = [sps_quants[sp] for sp in sps_1.keys()]
        vars_2 = [sps_quants[sp] for sp in sps_2.keys()]
        vals_1 = sps_1.values()
        vals_2 = sps_2.values()

        ratio_cons = bound_weighted_average_value_ratio(vars_1, vars_2, vals_1, vals_2, lb=lb, ub=ub)
        if return_constraint:
            return ratio_cons
        self.cons.append(ratio_cons)

    def synthesise_from(self, synthesis, return_constraint=False):
        return synthesis.fix_product(self.element_quantity_vars(), return_constraint)

    def bound_distance(self, other, calculator, return_constraint=False, *, ub=None, lb=None):
        """Constrain the distance between this composition and another given a calculator whose metric acts on normed element quantities.
        """
        elt_vars = {self.name: self.element_quantity_vars()}
        
        if isinstance(other, str):
            other = pg.Composition(other)

        if isinstance(other, pg.Composition):
            other = other.fractional_composition
            elt_vars.update({str(other): dict(other)})
        elif isinstance(other, TargetComposition):
            elt_vars.update({str(other): other.element_quantity_vars()})
        else:
            raise ValueError(f'Expected TargetComposition or pymatgen Composition object. Received type {type(other)}.')
        
        bound_dist_cons = calculator.bound_distance(elt_vars, lb, ub) # this relies on the returned constraints being ONLY those which represent the bound, not those that define the distance measure. 
        if return_constraint:
            return bound_dist_cons
        self.cons.append(bound_dist_cons)

    def property_predictor_category(self, model, n):
        vars = self.element_quantity_vars()
        vars = [(pg.Element(elt).Z, var) for elt, var in vars.items()]
        vars.sort(key=lambda x: x[0])
        vars = [v for k, v in vars]

        model.setup(vars)
        model.select_class(n)

    def format_solution(self, model, as_frac=False):
        out = {elt: model[elt_var] 
               for elt, elt_var in self.element_quantity_vars().items() 
               if model[elt_var].numerator_as_long() != 0}
        
        if as_frac:
            return {elt: Fraction(quant.numerator_as_long(), quant.denominator_as_long())
                    for elt, quant in out.items()}
        
        return {elt: round(float(quant.numerator_as_long()) / float(quant.denominator_as_long()), 3) 
                for elt, quant in out.items()}

    def __str__(self):
        return self.name
