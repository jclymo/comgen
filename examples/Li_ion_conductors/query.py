from comgen import SpeciesCollection, IonicComposition
from csv import DictReader 
import pymatgen.core as pg

reps_file = "LiIon_reps.csv"
ingredients_file = "Li_starting_materials.csv"
output_file = "results.txt"

distance = 1
num_results = 20

ingredients, excluded, comparisons = [], [], []

with open(ingredients_file) as f:
    for row in DictReader(f):
        comp = pg.Composition(row['norm_composition'])
        excluded.append(comp)
        if row['starting_material'] == 'Y':
            ingredients.append(comp)

with open(reps_file) as f:
    comparisons = [row['composition'] for row in DictReader(f)]

Li = {'Li'}
A = {'B', 'Al', 'Si', 'P'}
B = {'Mg', 'Zn', 'Ta', 'Sn', 'Ge', 'Ga'}
C = {'K', 'Ca', 'Sr', 'Y', 'Zr', 'Ba', 'La', 'Gd'}
D = {'N', 'O', 'F'}
E = {'S', 'Se', 'Te', 'Cl', 'Br', 'I'}

elts = Li|A|B|C|D|E
sps = SpeciesCollection.for_elements(elts)

query = IonicComposition(sps)

query.include_elements_quantity(Li, lb=0.2)
query.include_elements_quantity(A, lb=0.05)
query.include_elements_quantity(C, lb=0.05)
query.include_elements_quantity(D, lb=0.05)
query.include_elements_quantity(E, lb=0.05)

query.distinct_elements(ub=6)

query.elmd_far_from_all(comparisons, distance)

query.made_from(ingredients)

query.exclude(excluded)

query.total_atoms(lb=10, ub=15)

with open(output_file, 'w') as f_out:
    i = 0
    while i < num_results:
        res = query.get_next(as_frac=True)
        f_out.write(str(res)+'\n')
        i += 1

