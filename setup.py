from setuptools import setup

setup(
    name="comgen", 
    version="0.0.8",
    description="explore ionic composition space",
    packages=['comgen', 'comgen.constraintsystems', 'comgen.data'], 
    package_dir={'':'src'},
    package_data={'comgen.data': ['periodic_table.json', 'common_poly_ions.txt']},    
    install_requires=['pymatgen>=2022.5.26', 'z3-solver>=4.8.17.0']
)