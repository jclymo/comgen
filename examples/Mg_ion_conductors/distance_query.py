from comgen import SpeciesCollection, IonicComposition
from csv import DictReader 
import pymatgen.core as pg

li_conductors_file = "../data/LiIonDatabase.csv"
output_file = "results.txt"

distance = 2
num_results = 5

Mg = {'Mg'}
A = {'S', 'Se', 'Te', 'B', 'Al', 'Si', 'P', 'Zn', 'Ta', 'Sn', 'Ge', 'Ga', 'K', 'Ca', 'Sr', 'Y', 'Zr', 'Ba', 'La', 'Gd'}
B = {'N', 'O', 'F', 'Cl', 'Br', 'I'}
elts = Mg|A|B
sps = SpeciesCollection.for_elements(elts)

query = IonicComposition(sps, precision=0.01)
query.include_elements_quantity(Mg, lb=0.1)
query.distinct_elements(ub=6)
query.total_atoms(lb=10, ub=20)
with open(li_conductors_file) as f:
    comps = [row['composition'].strip('"') for row in DictReader(f) if float(row['target']) > float(5E-2)]
    query.elmd_close_to_one(comps, distance)

with open(output_file, 'w') as f_out:
    i = 0
    while i < num_results:
        res = query.get_next(as_frac=True)
        f_out.write(str(res)+'\n')
        i += 1




