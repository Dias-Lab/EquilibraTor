from setuptools import setup, find_packages

VERSION = '1.0.0'
DESCRIPTION = 'EquilibraTor'
setup(
    name="EquilibraTor",
    version=VERSION,
    author="José D. D. Cediel-Becerra",
    author_email="jcedielbecerra@ufl.edu",
    description=DESCRIPTION,
    packages=find_packages(),
    package_data={"": ["nvt_stage.mdp", "npt_stage.mdp", "ions.mdp", "minimization_stage.mdp", "production_stage.mdp"]},
    scripts=['bin/EquilibraTor'])
