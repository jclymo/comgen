from collections import defaultdict
from .data import data
import pymatgen.core as pg
import numpy as np

class PolyAtomicSpecies:
  """
    Mimics some features of pymatgen Species for compositions instead of elements.
    specifically, able to call species.oxi_state and species.X (for electronegativity)
  """
  def __init__(self, comp, oxi_state):
    if isinstance(comp, str):
      comp = pg.Composition(comp)

    self._composition = comp
    self._oxi_state = oxi_state 
    self._electronegativity = np.nan

  @property
  def oxi_state(self):
    return self._oxi_state

  @property
  def composition(self):
    return self._composition

  @property
  def X(self):
    return self._electronegativity
  
  def multiplier(self, element: str) -> int:
    """
      How many atoms of this element in the species? 
      e.g. for CO3, multiplier("C")=1, multiplier("O")=3.
      expects element as a string (not pymatgen element)
    """
    return int(dict(self.composition).get(element, 0))

  @property
  def elements(self):
    return set(self.composition.elements)
  
  def __str__(self) -> str:
    charge_str = '' if self.oxi_state in [1, -1] else str(abs(self.oxi_state))
    charge_str += '-' if self.oxi_state < 0 else '+'

    return "("+self.composition.to_pretty_string()+")"+charge_str

  def __repr__(self) -> str:
    charge_str = '' if self.oxi_state in [1, -1] else str(abs(self.oxi_state))
    charge_str += '-' if self.oxi_state < 0 else '+'

    return "PolyAtomicSpecies ("+self.composition.to_pretty_string()+")"+charge_str  
    
  def __eq__(self, other) -> bool:
    return self.composition == other.composition and self.oxi_state == other.oxi_state  

  def __hash__(self) -> int:
    return hash((self.composition, self.oxi_state))


class SpeciesCollection:
  def __init__(self, data: set):
    self._species = set()
    self.update(data) 

  def _set_empty_views(self):
    self._elements_view = defaultdict(set)
    self._compositions_view = defaultdict(set) 

  @staticmethod
  def _is_valid_collection_data(data):
    for item in data:
      assert isinstance(item, pg.Species) or isinstance(item, PolyAtomicSpecies)

  def ungrouped_view(self) -> set:
    return self._species

  # def group_by_pettifor_view(self):
  #   return {el.mendeleev_no: species for el, species in self.group_by_element_view().items()}

  def group_by_element_view(self):
    if self._elements_view:
      return self._elements_view
    
    for sp in self._species:
      if isinstance(sp, PolyAtomicSpecies):
        for el in sp.elements:
          self._elements_view[el].add(sp)
      else:
        self._elements_view[sp.element].add(sp)

    self._elements_view.default_factory = None # freeze default dict
    
    return self._elements_view

  def group_by_composition_view(self):
    if self._compositions_view:
      return self._compositions_view
    
    for sp in self._species:
      if isinstance(sp, PolyAtomicSpecies):
        for comp in sp.composition:
          self._compositions_view[comp].append(sp)
      else:
        self._compositions_view[sp.element].append(sp)
    
    self._elements_view.default_factory = None # freeze default dict

    return self._compositions_view

  def update(self, species):
    if isinstance(species, PolyAtomicSpecies) or isinstance(species, pg.Species):
      species = set(species) # to allow passing in a single item
    self._is_valid_collection_data(species) # error checking
    self._set_empty_views() # old derived data now invalid

    self._species.update(species)

  def difference(self, species):
    if isinstance(species, PolyAtomicSpecies) or isinstance(species, pg.Species):
      species = set(species) # to allow passing in a single item
    self._is_valid_collection_data(species) # error checking
    self._set_empty_views() # old derived data now invalid

    return SpeciesCollection(self._species.difference(species))

  def filter(self, elements: set):
    elements = {pg.Element(el) for el in elements}
    filtered_collection = self._species.copy()

    for item in self._species:
      if isinstance(item, pg.Species):
        if not item.element in elements:
          filtered_collection.remove(item)
      if isinstance(item, PolyAtomicSpecies):
        if not item.elements.issubset(elements):
          filtered_collection.remove(item)
    
    return SpeciesCollection(filtered_collection)

  def filter_mono_species(self):
    filtered_collection = set()
    
    for item in self._species:
      if isinstance(item, pg.Species):
        filtered_collection.add(item)

    return SpeciesCollection(filtered_collection)

  @classmethod
  def for_elements(
    cls, 
    elements=None, 
    permitted=data.PermittedSpecies.SH_FIX_BrClFINS, 
    include_poly=False):
    """
    For the given set of elements (or all elements), get collection of possible species.
    Specify whether only very common or all species are included.
    Specify whether polyatomic species (made of *only* these elements) are included.
    For fine control over included species can directly add / remove from the returned collection.
    """
    if elements:
      elements = {pg.Element(el) for el in elements}

    mono_species = data.get_permitted_mono_species(elements, permitted)
    species = {pg.Species(el, chg) for (el, chg) in mono_species}
    if include_poly:
      poly_species = data.get_poly_atomic_species(elements)
      species.update({PolyAtomicSpecies(comp, chg) for (comp, chg) in poly_species})

    return cls(species)

  def __len__(self):
    return len(self._species)

  def __str__(self):
    return f'{__class__.__name__} {self.ungrouped_view()}'





