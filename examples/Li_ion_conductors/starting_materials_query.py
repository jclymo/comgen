from comgen import SpeciesCollection, IonicComposition
from csv import DictReader 
import pymatgen.core as pg
from pathlib import Path

examples_dir = Path(__file__).resolve().parent.parent
data_dir = examples_dir / "data"
output_dir = examples_dir / "output"
output_file = output_dir / "li_starting_materials.txt"
reps_file = data_dir / "LiIon_reps.csv"
ingredients_file = data_dir / "Li_starting_materials.csv"

distance = 5
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

# remove unwanted species from default - for this query we only want to allow one charge per element
sps = sps.difference({
    pg.Species('Ta', 4), pg.Species('Ta', 3),
    pg.Species('Se', 4), pg.Species('Se', 6),
    pg.Species('Te', 4), pg.Species('Te', 6),
    pg.Species('Ge', 4),
    pg.Species('P', 3)})


query = IonicComposition(sps)

query.include_elements_quantity(Li, lb=0.2)

query.distinct_elements(ub=6)

query.elmd_far_from_all(comparisons, distance)

query.made_from(ingredients)

query.exclude(excluded)

query.total_atoms(lb=10, ub=15)

output_dir.mkdir(parents=True, exist_ok=True)
with open(output_file, "w", encoding="utf-8") as f_out:
    for _ in range(num_results):
        out = query.get_next(as_frac=True)
        if out is None:
            break
        res, _ = out
        f_out.write(str(res) + '\n')
