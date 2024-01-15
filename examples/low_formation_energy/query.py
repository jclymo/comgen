"""
minimal example querying a model which expects normed composition array as input.
This model predicts formation energy to be stable / metastable / unstable.
"""

from comgen import IonicComposition, SpeciesCollection
import onnx

model_path = 'ehull_1040_bn.onnx'
onnx_model = onnx.load(model_path)

sps = SpeciesCollection.for_elements()
query = IonicComposition(sps)
query.category_prediction(onnx_model, 0)

print(query.get_next())