from z3 import Real, Sum, And
import onnx
import numpy as np
from comgen.constraint_system.common import ReLU

class ONNX:
    def __init__(self, onnx_model, constraint_log):
        self.name = f'onnx{id(self)}'
        self.onnx_model = onnx_model
        self.initNames = [node.name for node in onnx_model.graph.initializer]
        self.inputNames = [inp.name for inp in onnx_model.graph.input if inp.name not in self.initNames]
        self.outputNames = [out.name for out in onnx_model.graph.output if out.name not in self.initNames]
        self._init_shapes()
        self.vars = {}
        self.cons = constraint_log

    def _init_shapes(self):
        self.shapes = {}
        for node in self.onnx_model.graph.input:
            dims = [d.dim_value for d in node.type.tensor_type.shape.dim]
            self.shapes[node.name] = dims
        self.shapes.update({node.name: node.dims for node in self.onnx_model.graph.initializer})

    def _get_graph_output_names(self):
        return [node.name for node in self.onnx_model.graph.output]

    def _get_node_input_names(self, node):
        return [input_name for input_name in node.input if input_name not in self.initNames]
    
    def _get_node_init_names(self, node):
        return [input_name for input_name in node.input if input_name in self.initNames]
    
    def _get_init_values(self, node_name):
        for init in self.onnx_model.graph.initializer:
            if init.name == node_name:
                data = np.frombuffer(init.raw_data, dtype=onnx.mapping.TENSOR_TYPE_MAP[1].np_dtype)
                data = data.reshape(self.shapes[node_name])
                return data

    def _get_attribute(self, node, att_name):
        for att in node.attribute:
            if att.name == att_name:
                return onnx.helper.get_attribute_value(att)  

    def _create_variables(self, node_name):
        assert node_name not in self.vars
        num_vars = self.shapes[node_name]
        self.vars[node_name] = np.array([Real(f'{self.name}_{node_name}_{i}') for i in range(np.prod(num_vars))]).reshape(num_vars)

    def _output_shape_after_broadcast(self, shape1, shape2):
        if shape1 == shape2:
            return shape1
        
        if len(shape2) < len(shape1): # want shape1 <= shape2
            shape1, shape2 = shape2, shape1
        diff = len(shape2)-len(shape1)
        shape1 = [1]*diff + shape1 # now shape1 == shape2
        return [max(s1, s1) for s1, s2 in zip(shape1, shape2)] # mimic broadcasting

    def make_gemm_equations(self, node):
        node_name = node.output[0]

        node_inputs = self._get_node_input_names(node)
        node_inits = self._get_node_init_names(node)

        assert len(node_inputs) == 1
        input = node_inputs[0]
        assert len(node_inits) == 2
        if len(self.shapes[node_inits[0]]) == 2 and len(self.shapes[node_inits[1]]) == 1:
            weight, bias = node_inits 
        elif len(self.shapes[node_inits[0]]) == 1 and len(self.shapes[node_inits[1]]) == 2:
            bias, weight = node_inits 
        else:
            assert False, f'cannot assign bias and weight inputs for Gemm layer {node.name}'

        in_shape = self.shapes[input]
        weight_shape = self.shapes[weight]
        bias_shape = self.shapes[bias]

        weight_data = self._get_init_values(weight)
        bias_data = self._get_init_values(bias)

        # TODO trans_a, alpha, beta
        tb = self._get_attribute(node, 'transB')
        if tb: 
            weight_shape = weight_shape[::-1]
            weight_data = weight_data.transpose()
        
        assert in_shape[-1] == weight_shape[0], f'{in_shape=}, {weight_shape=}' # in features
        assert weight_shape[1] == bias_shape[0] # out features
        dim_in, dim_out, dim_batch = weight_shape[0], weight_shape[1], np.prod(in_shape[:-1])

        self.shapes[node_name] = in_shape[:-1] + [weight_shape[1]]

        in_vars = self.vars[input]
        assert list(in_vars.shape) == in_shape, f'{in_vars.shape}, {in_shape}'
        in_vars = in_vars.reshape(-1, dim_in) # flatten batch dimensions, linear transformation on innermost dimension.
        
        self._create_variables(node_name)
        out_vars = self.vars[node_name].reshape(-1, dim_out)

        for i in range(dim_batch):
            for j in range(dim_out):
                gemm_con = Sum([in_vars[i][k] * weight_data[k][j] for k in range(dim_in)]) + bias_data[j] == out_vars[i][j]
                self.cons.append(gemm_con)

    def make_relu_equations(self, node):
        node_name = node.output[0]

        node_inputs = self._get_node_input_names(node)
        node_inits = self._get_node_init_names(node)
        assert len(node_inputs) == 1
        input = node_inputs[0]
        assert len(node_inits) == 0

        self.shapes[node_name] = self.shapes[input]
        dim_batch = np.prod(self.shapes[node_name])

        in_vars = self.vars[input]
        assert list(in_vars.shape) == self.shapes[node_name], f'{in_vars.shape}, {self.shapes[node_name]}'
        in_vars = in_vars.flatten() # flatten to act element-wise 

        self._create_variables(node_name)
        out_vars = self.vars[node_name].flatten()

        for i in range(dim_batch):
            relu_con = out_vars[i] == ReLU(in_vars[i])
            self.cons.append(relu_con)

    def make_add_equations(self, node):
        node_name = node.output[0]
        node_inputs = self._get_node_input_names(node)
        node_inits = self._get_node_init_names(node)
        assert len(node_inputs) + len(node_inits) == 2
        input1, input2 = node_inputs + node_inits
        
        self.shapes[node_name] = self._output_shape_after_broadcast(
            self.shapes[input1], 
            self.shapes[input2])

        self._create_variables(node_name)
        out_vars = self.vars[node_name].flatten()

        in_vars1 = self.vars[input1] if len(node_inputs) > 0 else self._get_init_values(input1)
        in_vars2 = self.vars[input2] if len(node_inputs) == 2 else self._get_init_values(input2)
        in_vars1 = np.broadcast_to(in_vars1, self.shapes[node_name]).flatten()
        in_vars2 = np.broadcast_to(in_vars2, self.shapes[node_name]).flatten()
        
        for i in range(len(out_vars)):
            add_con = out_vars[i] == in_vars1[i] + in_vars2[i]
            self.cons.append(add_con)

    def make_batchnorm_equations(self, node):       
        node_name = node.output[0]
        node_inputs = self._get_node_input_names(node)
        assert len(node_inputs) == 1
        input = node_inputs[0]
        node_inits = self._get_node_init_names(node)
        assert len(node_inits) == 4

        epsilon = self._get_attribute(node, "epsilon")
        weight = self._get_init_values(node.input[1]).flatten()
        bias = self._get_init_values(node.input[2]).flatten()
        mean = self._get_init_values(node.input[3]).flatten()
        variance = self._get_init_values(node.input[4]).flatten()

        in_vars = self.vars[input].flatten()

        assert len(weight) == len(in_vars) # can be False for CNNs but should be True for FF

        self.shapes[node_name] = self.shapes[input]
        self._create_variables(node_name)
        out_vars = self.vars[node_name].flatten()    

        for i in range(len(weight)):
            bn_con = out_vars[i] == (1 / np.sqrt(variance[i] + epsilon) * weight[i]) * in_vars[i] \
                - mean[i] / np.sqrt(variance[i] + epsilon) * weight[i] - bias[i]
            self.cons.append(bn_con)

    def add_equations(self, node):
        if node.op_type == 'Gemm':
            self.make_gemm_equations(node)
        elif node.op_type == 'Relu':
            self.make_relu_equations(node)
        elif node.op_type == 'Add':
            self.make_add_equations(node)
        elif node.op_type == 'BatchNormalization':
            self.make_batchnorm_equations(node)
        elif node.op_type == 'Sigmoid':
            assert node.output[0] in self.outputNames, 'Sigmoid only permitted at output.'
            # ignore -- constraints will be placed on previous layer
        elif node.op_type == 'Softmax':
            assert node.output[0] in self.outputNames, 'Softmax only permitted at output.'
            # ignore -- constraints will be placed on previous layer
        else:
            raise ValueError(f'Layer type {node.op_type} not yet implemented.')

    def get_node(self, node_name):
        node_list = [node for node in self.onnx_model.graph.node if node_name in node.output]
        assert len(node_list) == 1, f"{node_list=}, {node_name=}"
        return node_list[0]

    def process_node(self, node_name):
        if node_name in self.inputNames:
            return
        node = self.get_node(node_name)
        self.processed.append(node_name) 
        node_inputs = self._get_node_input_names(node)
        node_inits = self._get_node_init_names(node)
        for inp in node_inputs:
            if not inp in self.processed:
                self.process_node(inp)
        self.add_equations(node)       

    def setup(self, input_vars=None):
        if input_vars is not None:
            assert len(input_vars) == np.prod(self.shapes[self.inputNames[0]]) # assumes one input vector only
            node_name = self.inputNames[0]
            assert node_name not in self.vars
            num_vars = self.shapes[node_name]
            self.vars[node_name] = np.array(input_vars).reshape(num_vars)

        else:
            self._create_variables(self.inputNames[0])

        self.processed = []
        for out_node in self._get_graph_output_names():
            self.process_node(out_node)

    def select_class(self, n, return_constraint=False):
        model_output_names = self.outputNames[0]
        model_output_vars = self.vars[model_output_names][0]
        
        assert n < len(model_output_vars)
        
        cons = []
        target_var = model_output_vars[n]        
        for i, var in enumerate(model_output_vars):
            if i != n:
                cons.append(target_var > var)
        
        if return_constraint:
            return And(cons)
        self.cons.append(And(cons))

        # if final layer is sigmoid or softmax
        # constrain penultimate layer instead

    def fix_output(self, lb, ub, return_constraint=False):
        model_output_names = self.outputNames[0]
        model_output_vars = self.vars[model_output_names][0]

        assert len(lb) == len(model_output_vars)
        assert len(ub) == len(model_output_vars)

        cons = []
        for i, var in enumerate(model_output_vars):
            cons.append(var >= lb[i])
            cons.append(var <= ub[i])

        if return_constraint:
            return And(cons)
        self.cons.extend(cons)
        
        