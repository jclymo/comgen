from setuptools import setup

setup(
    name="comgen", 
    version="0.0.15",
    description="explore chemical compositions",
    packages=['comgen'],#, 'comgen.constraint_systems', 'comgen.constraint_systems.composition' 'comgen.data'], 
    # package_dir={'':'src'},
    # package_data={'comgen.data': ['periodic_table.json', 'common_poly_ions.txt']},    
    install_requires=['pymatgen>=2022.5.26', 'z3-solver>=4.8.17.0']
)