from z3 import sat, Solver

class Query:
    def __init__(self):
        self.constraints = []
        self.solutions = [] 

    def get_next(self):
        s = Solver() 
        for con in self.constraints:
            s.add(con)
        if s.check() != sat:
            return 
        model = s.model()
        self.solutions.append(model)
        return model
    