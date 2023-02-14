from enum import Enum
import pymatgen.core as pg
from pathlib import Path
import json
import csv
from ast import literal_eval

_pt_data_file = "periodic_table.json"
_poly_data_file = "common_poly_ions.txt"

with open(str(Path(__file__).absolute().parent / _pt_data_file)) as f:
  PT_DATA = json.load(f)
    
with open(str(Path(__file__).absolute().parent / _poly_data_file)) as f:
  ions_reader = csv.reader(f, delimiter='\t')
  POLY_DATA = {name: literal_eval(charges) for (name, charges) in ions_reader}

## TODO: add data files for other permitted species levels 

ELEMENTS = {pg.Element(el) for el in PT_DATA.keys() if pg.Element(el).Z <= 103}

class PermittedSpecies(Enum):
  FIXED = 1
  VARY_TRANSITION = 2
  SH_FIX_BrClFINS = 3
  SHANNON = 4
  ALL = 5

def get_permitted_mono_species(
  elements=None, 
  permitted=PermittedSpecies.SHANNON) -> set:
  """
  Get possible oxidation states for a given set of elements.
  """
  if elements is None:
    elements = ELEMENTS

  # if permitted == PermittedSpecies.FIXED:
  #   species = {(el, chg) for el in elements for chg in FIXED_DATA[el] }
  
  if permitted == PermittedSpecies.SHANNON:
    el_radii = {el: PT_DATA[el.symbol].get('Shannon radii', {"0":{}}) for el in elements}
    species = {(el, int(ch)) for el, item in el_radii.items() for ch in item.keys()}
  
  if permitted == PermittedSpecies.SH_FIX_BrClFINS:
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

def get_poly_atomic_species(elements=None) -> set:
  if elements is None:
    elements = ELEMENTS
  
  species = set()

  for c, chgs in POLY_DATA.items():
    comp = pg.Composition(c)
    if set(comp.elements).issubset(set(elements)):
      species.update({(comp, chg) for chg in chgs})

  return species
