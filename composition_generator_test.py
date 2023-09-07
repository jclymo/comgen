from comgen import SpeciesCollection, IonicCompositionGenerator
import comgen as cc
'''
elts = {'Mg', 
'B', 'Al', 'Si', 'P', 'K', 'Ca', 'Zn', 'Sr', 'Y', 'Zr', 'Sn', 'Ba', 'Ta', 'La', 
'N', 'O', 'S', 'F', 'Cl', 'Br', 'I'}

# get all valid species for the elements we're interested in
sps = SpeciesCollection.for_elements(elts)

# create the generator
gen = IonicCompositionGenerator(sps)

# pick exactly four elements from the following sets
gen.distinct_elements(4)

gen.include_element_from({'Mg'})
gen.include_element_from({'B', 'Al', 'Si', 'P', 'K', 'Ca', 'Zn', 
                        'Sr', 'Y', 'Zr', 'Sn', 'Ba', 'Ta', 'La'}, lb=1, ub=2)
gen.include_element_from({'N', 'O', 'S', 'F', 'Cl', 'Br', 'I'}, lb=1, ub=2)
# print(gen.get_next(2))


# Mg at least 0.1 of total
gen.fix_elements_quantity('Mg', lb=0.1)
gen.fix_elements_quantity({'B', 'Al', 'Si', 'P', 'K', 'Ca', 'Zn', 
                        'Sr', 'Y', 'Zr', 'Sn', 'Ba', 'Ta', 'La'}, lb=0.2)

print(gen.get_next())

# add comparison compositions
gen.emd_comparison_compositions(
    [
        'Li3.4Si0.4P0.6S4', 
        'Li6.4La3Zr1.4Ta0.6O12', 
        # 'Li7.06La3Zr1.94Y0.06O12', 
        # 'Li6.4Al0.19La3Zr2O11.8', 
        # 'La0.57Li0.29TiO3', 
        # 'Li6.5La3Zr1.75Te0.25O12'
    ], 

    ub=3)
print(gen.get_next())

# force smaller unit cell
gen.total_atoms(lb=10, ub=20)
print(gen.get_next())

# print(gen.get_next(5))

print('\n'.join(gen.get_constraints_summary()))


'''

import pymatgen.core as pg
from comgen.data import data
import unittest

class IonComGenTest(unittest.TestCase):
    elts = {'Li', 
    'B', 'Al', 'Si', 'P',
    'Mg', 'Zn', 'Ta', 'Sn', 'Ge', 'Ga',
    'K', 'Ca', 'Sr', 'Y', 'Zr', 'Ba', 'La', 'Gd', 
    'N', 'O', 'F',
    'S', 'Se', 'Te', 'Cl', 'Br', 'I'}
    def __init__(self, verbose=False):
        super().__init__()

        # self.sps = cc.SpeciesCollection.for_elements(self.elts, permitted=data.PermittedSpecies.SHANNON, include_poly=False)
        self.reset_sps()
        self.verbose = verbose

    def reset_sps(self):
        sps = SpeciesCollection.for_elements(self.elts, include_poly=False)
        self.sps = sps.difference({
            pg.Species('Ta', 4), pg.Species('Ta', 3), 
            pg.Species('Se', 4), pg.Species('Se', 6),
            pg.Species('Te', 4), pg.Species('Te', 6),
            pg.Species('Ge', 4),
            pg.Species('P', 3),
            })

    def charge_balance(self):
        def fail_eneg():
            if self.verbose: print('testing charge balance: fix quantities Mg0.5 O0.5 with Mg-2 and O+2 (cannot balance charge and respect electronegativity) (UNSAT)')
            self.sps.update({pg.Species('Mg', -2), pg.Species('O', 2)})
            self.sps = self.sps.difference({pg.Species('Mg', 2), pg.Species('O', -2)})
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.fix_elements_quantity('Mg', 0.5)
            gen.fix_elements_quantity('O', 0.5)
            return gen.get_next()

        def fail_balance():
            if self.verbose: print('testing charge balance: fix quantities Al0.5 O0.5 with standard charges (UNSAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.fix_elements_quantity('Al', 0.5)
            gen.fix_elements_quantity('O', 0.5)
            return gen.get_next()
        
        def fail_precision():
            if self.verbose: print('testing charge balance: fix quantities Cl0.286 O0.714 with standard charges (UNSAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            # this doesn't work because not sufficiently accurate. 
            gen.fix_elements_quantity('Cl', 0.286) #5
            gen.fix_elements_quantity('O', 0.714) #-2            
            return gen.get_next()

        def success():
            if self.verbose: print('testing charge balance: fix quantities Ta 2/7 O 5/7 with standard charges (SAT)') 
            gen = cc.IonicCompositionGenerator(self.sps)
            from z3 import Q # rational number
            gen.fix_elements_quantity('Ta', Q(2,7)) #5
            gen.fix_elements_quantity('O', Q(5,7)) #-2
            return gen.get_next()

        self.assertNotEqual(success(), [])
        self.assertEqual(fail_balance(), [])
        self.assertEqual(fail_precision(), [])
        self.assertEqual(fail_eneg(), []) 
        self.reset_sps() # because the previous test changed the permitted charges. 

    def distinct_elements(self):
        def success():
            if self.verbose: print('testing distinct_elements: fix quantities Mg0.25O0.5Zn0.125Ca0.125 and require 4 distinct elements (SAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.distinct_elements(4)
            gen.fix_elements_quantity('Mg', 0.25)
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Zn', 0.125)
            gen.fix_elements_quantity('Ca', 0.125)
            return gen.get_next()

        def fail():
            if self.verbose: print('testing distinct_elements: fix quantities Mg0.25O0.5Zn0.125Ca0.125 and require 5 distinct elements (UNSAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.distinct_elements(4)    
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Mg', 0.125)
            gen.fix_elements_quantity('Zn', 0.125)
            gen.fix_elements_quantity('Ca', 0.125)
            gen.fix_elements_quantity('Sr', 0.125)
            return gen.get_next()
        
        self.assertNotEqual(success(), [])   
        self.assertEqual(fail(), [])
    
    def select_from_set(self):
        def success():
            if self.verbose: print('testing include elt from set: fix quantities Mg0.25O0.5Zn0.125Ca0.125 and require include from {"Mg", "Ta"} (SAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.include_element_from({'Mg', 'Ta'})
            gen.fix_elements_quantity('Mg', 0.25)
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Zn', 0.125)
            gen.fix_elements_quantity('Ca', 0.125)
            return gen.get_next()
        
        def fail():
            if self.verbose: print('testing include elt from set: fix quantities Sr0.25O0.5Zn0.125Ca0.125 and require include from {"Mg", "Ta"} (UNSAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.include_element_from({'Mg', 'Ta'})
            gen.fix_elements_quantity('Sr', 0.25)
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Zn', 0.125)
            gen.fix_elements_quantity('Ca', 0.125)
            return gen.get_next()            
        
        self.assertNotEqual(success(), [])   
        self.assertEqual(fail(), [])

    def fix_elements_quantity(self):
        def success():
            if self.verbose: print('testing fix_elements_quantity: fix quantities Mg0.25O0.5Zn0.125Ca0.125 and require Mg quantity between 0.1 and 1 (SAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.fix_elements_quantity('Mg', lb=0.1, ub=1)
            gen.fix_elements_quantity('Mg', 0.25)
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Zn', 0.125)
            gen.fix_elements_quantity('Ca', 0.125)
            return gen.get_next()
        
        def fail():
            if self.verbose: print('testing fix_elements_quantity: fix quantities Mg0.25O0.5Zn0.125Ca0.125 and require Mg quantity between 0.3 and 1 (UNSAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.fix_elements_quantity('Mg', lb=0.3, ub=1)
            gen.fix_elements_quantity('Mg', 0.25)
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Zn', 0.125)
            gen.fix_elements_quantity('Ca', 0.125)
            return gen.get_next()

        self.assertNotEqual(success(), [])   
        self.assertEqual(fail(), [])

    def emd(self):
        def success():
            if self.verbose: print('testing emd: fix quantities Mg0.17O0.5Zn0.25Ca0.08 and require MgOCaZn between 16.65 and 16.7 (SAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.emd_comparison_compositions(['MgOCaZn'], lb=16.65, ub=16.7) # should be 16.69
            gen.fix_elements_quantity('Mg', 0.17)
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Zn', 0.25)
            gen.fix_elements_quantity('Ca', 0.08)
            return gen.get_next()
        
        def fail():
            if self.verbose: print('testing emd: fix quantities Mg0.17O0.5Zn0.25Ca0.08 and require MgOCaZn between 0 and 16.6 (UNSAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.emd_comparison_compositions(['MgOCaZn'], lb=0, ub=16.6) # should be 16.69
            gen.fix_elements_quantity('Mg', 0.17)
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Zn', 0.25)
            gen.fix_elements_quantity('Ca', 0.08)
            return gen.get_next()            

        self.assertNotEqual(success(), [])
        self.assertEqual(fail(), [])
        # TODO test with multiple compositions


    def exclude_solution(self):
        def success():
            if self.verbose: print('testing exclude_solution: fix quantities Mg0.17 O0.5 Ca0.08 and find two solutions (SAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.fix_elements_quantity('Mg', 0.17)
            gen.fix_elements_quantity('O', 0.5)
            # gen.fix_element_quantity('Zn', 0.25)
            gen.fix_elements_quantity('Ca', 0.08)
            gen.get_next()
            return gen.get_next()

        def fail():
            if self.verbose: print('testing exclude_solution: fix quantities Mg0.17 O0.5 Zn0.25 Ca0.08 and find two solutions (UNSAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.fix_elements_quantity('Mg', 0.17)
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Zn', 0.25)
            gen.fix_elements_quantity('Ca', 0.08)
            gen.get_next()
            return gen.get_next()

        self.assertNotEqual(success(), [])
        self.assertEqual(fail(), [])

    def total_atoms(self):
        def success():
            if self.verbose: print('testing max total atoms: fix quantities Mg0.25 O0.5 Zn0.125 Ca0.125 and require at most 8 atoms (SAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.max_total_atoms(8)
            gen.fix_elements_quantity('Mg', 0.25)
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Zn', 0.125)
            gen.fix_elements_quantity('Ca', 0.125)
            return gen.get_next()

        def fail():
            if self.verbose: print('testing max total atoms: fix quantities Mg0.25 O0.5 Zn0.125 Ca0.125 and require at most 7 atoms (UNSAT)')
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.max_total_atoms(7)
            gen.fix_elements_quantity('Mg', 0.25)
            gen.fix_elements_quantity('O', 0.5)
            gen.fix_elements_quantity('Zn', 0.125)
            gen.fix_elements_quantity('Ca', 0.125)
            return gen.get_next()
        
        self.assertNotEqual(success(), [])   
        self.assertEqual(fail(), [])    

    def construct_from(self):
        def success():
            if self.verbose: print("testing construct from compositions: fix quantities Li5 Cl1 O5 Al2 and input comps ['Cl2Li2', 'Al4O6', 'Li2O1'] (SAT)")
            from z3 import Q
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.fix_elements_quantity('Li', Q(5, 13))
            gen.fix_elements_quantity('Cl', Q(1, 13))
            gen.fix_elements_quantity('O', Q(5, 13))
            gen.fix_elements_quantity('Al', Q(2, 13))
            gen.construct_from(['Cl2Li2', 'Al4O6', 'Li2O1'])
            return gen.get_next()
        
        def fail():
            if self.verbose: print("testing construct from compositions: fix quantities Li5 Cl1 O5 Al2 and input comps ['Cl2Li2', 'Al1O1Cl1', 'Li2O1'] (UNSAT)")
            from z3 import Q
            gen = cc.IonicCompositionGenerator(self.sps)
            gen.fix_elements_quantity('Li', Q(5, 13))
            gen.fix_elements_quantity('Cl', Q(1, 13))
            gen.fix_elements_quantity('O', Q(5, 13))
            gen.fix_elements_quantity('Al', Q(2, 13))
            gen.construct_from(['Cl2Li2', 'Al1O1Cl1', 'Li2O1'])
            return gen.get_next()

        self.assertNotEqual(success(), [])
        self.assertEqual(fail(), [])

    def run(self):
        if self.verbose: print('starting tests.')
        self.charge_balance()
        self.distinct_elements()
        self.fix_elements_quantity()
        self.select_from_set()
        self.emd()
        self.exclude_solution()
        self.total_atoms()
        self.construct_from()
        print('all tests passed.')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    "--verbose",
    action="store_true",
    help="Output description of tests.",
)
args = parser.parse_args()

test = IonComGenTest(args.verbose)
test.run()
