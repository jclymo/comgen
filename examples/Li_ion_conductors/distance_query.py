from comgen import SpeciesCollection, IonicComposition
from csv import DictReader 
import pymatgen.core as pg
from pathlib import Path

output_file = "results.txt"

data_dir = Path(__file__).parent.parent / 'data'
reps_file = data_dir / 'LiIon_reps.csv'

distance = 5
num_results = 20

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
sps = sps.difference({
    pg.Species('Ta', 4), pg.Species('Ta', 3),
    pg.Species('Se', 4), pg.Species('Se', 6),
    pg.Species('Te', 4), pg.Species('Te', 6),
    pg.Species('Ge', 4),
    pg.Species('P', 3)})

query = IonicComposition(sps)

query.include_elements_quantity(Li, lb=0.2)
query.include_elements_quantity(A, lb=0.05)
query.include_elements_quantity(B, lb=0.05)
query.include_elements_quantity(C, lb=0.05)
query.include_elements_quantity(D, lb=0.05)
query.include_elements_quantity(E, lb=0.05)

query.distinct_elements(ub=6)

query.elmd_far_from_all(comparisons, distance)

query.total_atoms(lb=10, ub=15)

with open(output_file, 'w') as f_out:
    i = 0
    while i < num_results:
        res, model = query.get_next(as_frac=True)
        f_out.write(str(res)+'\n')
        i += 1
