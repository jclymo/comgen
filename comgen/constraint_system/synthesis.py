from comgen.constraint_system.common import zero_weighted_sum, weighted_sum
from z3 import Real, And

class Synthesis:
    def __init__(self, ingredient_compositions, constraint_log, return_vars):
        self.ingredient_compositions = ingredient_compositions
        self.name = f"Synthesis{id(self)}"
        self.cons = constraint_log
        self.return_vars = return_vars
        self._ingredient_quantity_variable_collection = {}
        self._setup()

    def _ingredient_quantity_vars(self, comp):
        return self._ingredient_quantity_variable_collection.get(str(comp))
    
    def _new_ingredient_quantity_var(self, comp):
        var = Real(f'{self.name}_{str(comp)}_ingredientquantity')
        self._ingredient_quantity_variable_collection[str(comp)] = var
        self.return_vars.append(var)
        return var

    def _setup(self):
        for comp in self.ingredient_compositions:
            self._new_ingredient_quantity_var(comp)

    def fix_product(self, element_quantities, return_constraint=False):
        """
        Given some input starting materials (ingredients) require that the final composition can be made from some combination of these ingredients.
        i.e. weighted (by amount of each ingredient) sum of ingredient compositions equals final composition (specified by provided element quantities)
        """
        elements = {str(el) for comp in self.ingredient_compositions for el in comp}
        elements.update(set(element_quantities.keys()))
        ing_vars = [self._ingredient_quantity_vars(comp) for comp in self.ingredient_compositions]
        
        cons = []
        for el in elements:
            vars = ing_vars + [element_quantities.get(el, 0)]
            weights = [comp.get(el, 0) for comp in self.ingredient_compositions]
            weights.append(-1)
            cons.append(zero_weighted_sum(vars, weights))

        if return_constraint:
            return And(cons)
        self.cons.append(And(cons))

    def bound_cost(self, ingredient_costs, lb, ub):
        assert all([k in self.ingredient_compositions for k in ingredient_costs.keys()])
        vars = self._ingredient_quantity_vars()
        weights = [ingredient_costs[comp] for comp in self.ingredient_compositions]
        self.cons.append(weighted_sum(vars, weights, lb=lb, ub=ub))