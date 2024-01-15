"""
Generate compositions which may have garnet-like structure.
Restrict the relative size of +2 and +3 cations.
Use SiO4 as the only anion. 
"""

from comgen import IonicComposition, SpeciesCollection, PolyAtomicSpecies

def get_permitted_ions():
    sps = SpeciesCollection.for_elements()
    sps = sps.having_charge([2,3])
    sps.update(PolyAtomicSpecies("SiO4", -4))
    return sps
    
sps = SpeciesCollection.for_elements() # all mono species
sps1 = sps.having_charge(2)
sps2 = sps.having_charge(3)

sps = get_permitted_ions()

query = IonicComposition(sps)

query.distinct_elements(4)
query.include_elements_quantity({'O'}, 12/20)
query.radius_ratio(sps.having_charge(2), sps.having_charge(3), cn1='VIII', cn2='VI', lb=1.5, ub=1.9)
print(query.get_next())
