import pymatgen.core as pg

PETTIFOR_KEYS = tuple(range(0, 103))

MOD_PETTI = {"D": 102, "T": 102, "H": 102, "He": 0, "Li": 11, 
"Be": 76, "B": 85, "C": 86, "N": 87, "O": 96, "F": 101, 
"Ne": 1, "Na": 10, "Mg": 72, "Al": 77, "Si": 84, "P": 88, 
"S": 95, "Cl": 100, "Ar": 2, "K": 9, "Ca": 15, "Sc": 47, 
"Ti": 50, "V": 53, "Cr": 54, "Mn": 71, "Fe": 70, "Co": 69, 
"Ni": 68, "Cu": 67, "Zn": 73, "Ga": 78, "Ge": 83, "As": 89, 
"Se": 94, "Br": 99, "Kr": 3, "Rb": 8, "Sr": 14, "Y": 20, "Zr": 48, 
"Nb": 52, "Mo": 55, "Tc": 58, "Ru": 60, "Rh": 62, "Pd": 64, "Ag": 66, 
"Cd": 74, "In": 79, "Sn": 82, "Sb": 90, "Te": 93, "I": 98, "Xe": 4, 
"Cs": 7, "Ba": 13, "La": 31, "Ce": 30, "Pr": 29, "Nd": 28, "Pm": 27, 
"Sm": 26, "Eu": 16, "Gd": 25, "Tb": 24, "Dy": 23, "Ho": 22, "Er": 21, 
"Tm": 19, "Yb": 17, "Lu": 18, "Hf": 49, "Ta": 51, "W": 56, "Re": 57, 
"Os": 59, "Ir": 61, "Pt": 63, "Au": 65, "Hg": 75, "Tl": 80, "Pb": 81, 
"Bi": 91, "Po": 92, "At": 97, "Rn": 5, "Fr": 6, "Ra": 12, "Ac": 32, 
"Th": 33, "Pa": 34, "U": 35, "Np": 36, "Pu": 37, "Am": 38, "Cm": 39, 
"Bk": 40, "Cf": 41, "Es": 42, "Fm": 43, "Md": 44, "No": 45, "Lr": 46, 
"Rf": 0, "Db": 0, "Sg": 0, "Bh": 0, "Hs": 0, "Mt": 0, "Ds": 0, "Rg": 0, 
"Cn": 0, "Nh": 0, "Fl": 0, "Mc": 0, "Lv": 0, "Ts": 0, "Og": 0, "Uue": 0}

def element_to_pettifor(elt):
    if isinstance(elt, str):
        elt = pg.Element(elt)
    assert isinstance(elt, pg.Element), elt
    return MOD_PETTI[elt.symbol]

def get_radii(sps, cn=None):
    radii = {}
    if cn is not None:
        for sp in sps.ungrouped_view():
            try:
                radii[str(sp)] = sp.get_shannon_radius(cn=cn, spin='High Spin', radius_type='crystal')
            except KeyError:
                pass
    else:
        radii = {str(sp): sp.ionic_radius for sp in sps.ungrouped_view()}
    return radii


