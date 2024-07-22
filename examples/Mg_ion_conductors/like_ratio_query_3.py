from comgen import SpeciesCollection, IonicComposition
from csv import DictReader 
import pymatgen.core as pg
from fractions import Fraction

num_results = 3
output_file = 'results.txt'

def like_Li6PS5Cl(denom):
    Mg = {'Mg'}
    A = {'S', 'Se', 'Te', 'B', 'Al', 'Si', 'P', 'Zn', 
        'Ta', 'Sn', 'Ge', 'Ga', 'K', 'Ca', 'Sr', 'Y', 
        'Zr', 'Ba', 'La', 'Gd', 'Mn', 
        'N', 'O', 'F', 'Cl', 'Br', 'I'}

    sps = SpeciesCollection.for_elements(A)
    mg_sps = SpeciesCollection.for_elements(Mg)

    mg_sps = mg_sps.having_charge(2)
    sps1 = sps.having_charge({1,2,3,4,5,6,7,8,9})
    sps2 = sps.having_charge({-1,-2,-3,-4,-5,-6,-7,-8,-9})
    sps.update(mg_sps)

    query = IonicComposition(sps, precision=0.01)

    query.ion_pair_radius_ratio(mg_sps, sps1, lb=1.55, ub=1.85)
    query.ion_pair_radius_ratio(mg_sps, sps2, lb=0.45, ub=0.55)

    query.distinct_elements(lb=3, ub=6)

    query.total_atoms(denom)

    query.include_elements_quantity(Mg, lb=Fraction(3, denom), ub=Fraction(6, denom))

    query.include_species_quantity(sps1, Fraction(1, denom))
    query.include_species_quantity(sps2, Fraction(6, denom)) 

    with open(output_file, 'a') as f_out:
        i = 0
        while i < num_results:
            res, model = query.get_next(as_frac=True)
            f_out.write(str(res)+'\n')
            i += 1

for n in range(10, 14):
    like_Li6PS5Cl(n)
