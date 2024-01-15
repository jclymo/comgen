# import pymatgen.core as pg
from comgen.constraint_system import ONNX
from comgen.query import Query
import unittest
import onnx
from z3 import Solver, sat
from numpy import argmax

class NNQuery(Query):
    def __init__(self, onnx_model):
        super().__init__()
        self.nn = ONNX(onnx_model, self.constraints)
        self.nn.setup()

    def fix_input(self, input):
        input_vars = self.nn.inputNames[0]
        assert len(self.nn.vars[input_vars][0]) == len(input)
        for var, val in zip(self.nn.vars[input_vars][0], input):
            self.constraints.append(var == val)      

    def solve_for_category(self, n):
        self.nn.select_class(n)
        model = self.get_next()
        return model is not None

    def solve_for_output(self, lb, ub):
        self.nn.fix_output(lb, ub)
        model = self.get_next()
        return model is not None

    def calculate_output(self, input):
        input_vars = self.nn.inputNames[0]
        assert len(self.nn.vars[input_vars][0]) == len(input)
        for var, val in zip(self.nn.vars[input_vars][0], input):
            self.constraints.append(var == val)

        model = self.get_next()

        model_output_names = self.nn.outputNames[0]
        model_output_vars = self.nn.vars[model_output_names][0]
        model_output_values = []
        for var in model_output_vars:
            val = model[var]
            val = round(float(val.numerator_as_long()) / float(val.denominator_as_long()), 4)
            model_output_values.append(val)   

        return model_output_values     


class ONNXTest(unittest.TestCase):
    def __init__(self, onnx_model):
        super().__init__()
        self.onnx_model = onnx_model

    def forward_pass(self, input, expected_output):
        query = NNQuery(self.onnx_model)
        output = query.calculate_output(input)
        self.assertEqual(output, expected_output)
    
    def select_class(self, input, n):
        query = NNQuery(self.onnx_model)
        query.fix_input(input)
        self.assertTrue(query.solve_for_category(n))

    def select_wrong_class(self, input, n):
        query = NNQuery(self.onnx_model)
        query.fix_input(input)
        self.assertFalse(query.solve_for_category(n))

    def bound_output(self, input, lb, ub):
        query = NNQuery(self.onnx_model)
        query.fix_input(input)    
        self.assertTrue(query.solve_for_output(lb, ub))    

def do_tests(network_path=None, test_input=None, expected_output=None):
    network_path = network_path or 'test_model.onnx'
    test_input = test_input or [0]*4
    expected_output = expected_output or [-0.2051,  0.2430,  0.0456]

    with open(network_path, "rb") as f:
        onnx_model = onnx.load(f)
    tester = ONNXTest(onnx_model)
    
    tester.forward_pass(test_input, expected_output)

    if len(expected_output) > 1:
        n = argmax(expected_output)
        m = 0 if n != 0 else len(expected_output)-1
        tester.select_class(test_input, n)
        tester.select_wrong_class(test_input, m)

    lb = [out - 0.1 for out in expected_output]
    ub = [out + 0.1 for out in expected_output]
    tester.bound_output(test_input, lb, ub)

if __name__ == "__main__":
    do_tests()







