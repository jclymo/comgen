from z3 import sat, Solver, Q
from fractions import Fraction

class Query:
    def __init__(self):
        self.constraints = []
        self.solutions = [] 

    def frac_to_rational(self, val):
        if isinstance(val, Fraction):
            return Q(val.numerator, val.denominator)
        return val

    def get_next(self):
        s = Solver() 
        for con in self.constraints:
            s.add(con)
        if s.check() != sat:
            return 
        model = s.model()
        self.solutions.append(model)
        return model
    