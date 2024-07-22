from comgen import SpeciesCollection, IonicComposition
from csv import DictReader 
from pathlib import Path

output_file = "results.txt"

data_dir = Path(__file__).parent.parent / 'data'
li_conductors_file = data_dir / 'LiIonDatabase.csv'

distance = 3
num_results = 5

Mg = {'Mg'}
A = {'S', 'Se', 'Te', 'B', 'Al', 'Si', 'P', 'Zn', 
     'Ta', 'Sn', 'Ge', 'Ga', 'K', 'Ca', 'Sr', 'Y', 
     'Zr', 'Ba', 'La', 'Gd',
     'N', 'O', 'F', 'Cl', 'Br', 'I'}
elts = Mg|A
sps = SpeciesCollection.for_elements(elts)

query = IonicComposition(sps, precision=0.01)
query.include_elements_quantity(Mg, lb=0.1)
query.distinct_elements(ub=6)
query.total_atoms(lb=10, ub=20)
comps = []
with open(li_conductors_file) as f:
    for row in DictReader(f):
        if float(row['target']) >= float(1E-3):
            if float(row['temperature']) >= 15 and float(row['temperature']) <= 35:
                comps.append(row['composition'].strip('"'))
    query.elmd_close_to_one(comps, distance)

with open(output_file, 'w') as f_out:
    i = 0
    while i < num_results:
        res, model = query.get_next(as_frac=True)
        f_out.write(str(res)+'\n')
        i += 1
