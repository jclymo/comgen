# import pymatgen.core as pg
from comgen.constraint_system import EMD
from comgen.query import Query
import unittest
from z3 import Solver, sat, Real

class EMDQuery(Query):
    def __init__(self):
        super().__init__()
        self.mapping_func = lambda x: x
        self.metric_ids = [0,1,2]
        self.vars = {1: {i: Real(f'{i}') for i in range(len(self.metric_ids))}}
        self.calculator = EMD(self.mapping_func, self.metric_ids, self.constraints)

    def setup(self, obj1, obj2):
        self.calculator._setup_distance_calculation({1:obj1, 2:obj2})

    def fix_input(self, vars, vals):
        for i in range(len(self.metric_ids)):
            self.constraints.append(vars[i] ==  vals[i])

    def calculate_output(self):
        model = self.get_next()
        output_name = self.calculator._distance_var(1,2)
        return model[output_name]        

class EMDTest(unittest.TestCase):
    def __init__(self):
        super().__init__()

    def forward(self, input1, input2, expected_output):
        query = EMDQuery()
        m_ids = list(range(len(input1)))
        vars = {i: Real(f'var_{i}') for i in m_ids}
        values = {i: input1[i] for i in m_ids}
        comparison = {i: input2[i] for i in m_ids}
        query.setup(vars, comparison)
        query.fix_input(vars, values)
        output = query.calculate_output()
        self.assertEqual(output, expected_output)   

def do_tests():
    tester = EMDTest()
    tester.forward([1,3,2], [3,2,1], 3)

if __name__ == "__main__":
    do_tests()







