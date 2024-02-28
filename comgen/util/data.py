from enum import Enum
import pymatgen.core as pg
import json
import csv
from ast import literal_eval
import pkg_resources

_pt_data_file = "periodic_table.json"
_poly_data_file = "common_poly_ions.txt"

with open(pkg_resources.resource_filename('comgen.util.data_files', _pt_data_file)) as f:
    PT_DATA = json.load(f)

with open(pkg_resources.resource_filename('comgen.util.data_files', _poly_data_file)) as f:
    ions_reader = csv.reader(f, delimiter='\t')
    POLY_DATA = {name: literal_eval(charges) for (name, charges) in ions_reader}

ELEMENTS = {pg.Element(el) for el in PT_DATA.keys() if pg.Element(el).Z <= 103}

class PossibleSpecies(Enum):
    FIXED = 1
    VARY_TRANSITION = 2
    SH_FIX_HALOGEN = 3 # Br- Cl- F- I- N-3 S-2  
    SHANNON = 4
    ALL = 5

def mono_atomic_species(
    elements=None, 
    permitted=PossibleSpecies.SHANNON) -> set:
    """
    Get possible oxidation states for a given set of elements.
    """
    if elements is None:
        elements = ELEMENTS

    if permitted == PossibleSpecies.SHANNON:
        el_radii = {el: PT_DATA[el.symbol].get('Shannon radii', {"0":{}}) for el in elements}
        species = {(el, int(ch)) for el, item in el_radii.items() for ch in item.keys()}

    if permitted == PossibleSpecies.SH_FIX_HALOGEN:
        el_radii = {el: PT_DATA[el.symbol].get('Shannon radii', {"0":{}}) for el in elements}
        species = set()
        for el, item in el_radii.items():
            if el.symbol in ['Br', 'Cl', 'F', 'I']:
                species.add((el, -1))
            elif el.symbol == 'N':
                species.add((el, -3))
            elif el.symbol == 'S':
                species.add(('S', -2))
            else:
                for ch in item.keys():
                    species.add((el, int(ch)))

    return species 

def poly_atomic_species(elements=None) -> set:
    if elements is None:
        elements = ELEMENTS

    species = set()

    for c, chgs in POLY_DATA.items():
        comp = pg.Composition(c)
        if set(comp.elements).issubset(set(elements)):
            species.update({(comp, chg) for chg in chgs})

    return species
