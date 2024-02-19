from comgen import SpeciesCollection, IonicComposition
from csv import DictReader 
import pymatgen.core as pg
from fractions import Fraction

num_results = 10

def like_Li6PS5Cl(output_file="results.txt", mg_lb=None, mg_ub=None):
    Mg = {'Mg'}
    A = {'S', 'Se', 'Te', 'B', 'Al', 'Si', 'P', 'Zn', 'Ta', 'Sn', 'Ge', 'Ga', 'K', 'Ca', 'Sr', 'Y', 'Zr', 'Ba', 'La', 'Gd'}
    B = { 'O', 'F', 'Cl', 'Br', 'I', 'N'}

    elts = A|B
    sps = SpeciesCollection.for_elements(elts)
    mg_sps = SpeciesCollection.for_elements(Mg)

    mg_sps = mg_sps.having_charge(2)
    sps1 = sps.having_charge({1,2,3,4,5,6,7,8,9}).difference(pg.Species('Mg', 2))
    sps2 = sps.having_charge({-1,-2,-3,-4,-5,-6,-7,-8,-9})
    sps.update(mg_sps)

    query = IonicComposition(sps, precision=0.01)

    query.ion_pair_radius_ratio(mg_sps, sps1, lb=1.55, ub=1.85)
    query.ion_pair_radius_ratio(mg_sps, sps2, lb=0.45, ub=0.55)

    query.distinct_elements(lb=3, ub=6)

    query.total_atoms(13)

    if mg_lb or mg_ub:
        query.include_elements_quantity(Mg, lb=Fraction(mg_lb, 13), ub=Fraction(mg_ub, 13))
    else:
        query.include_elements_quantity(Mg, Fraction(6,13))

    query.include_species_quantity(sps1, lb=Fraction(1,13))
    query.include_species_quantity(sps2, lb=Fraction(6,13)) 

    with open(output_file, 'a') as f_out:
        i = 0
        while i < num_results:
            res = query.get_next(as_frac=True)
            f_out.write(str(res)+'\n')
            i += 1


like_Li6PS5Cl('res1.txt')
like_Li6PS5Cl('res2.txt', 4, 5)
