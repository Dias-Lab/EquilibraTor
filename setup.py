from setuptools import setup, find_packages

VERSION = '0.5.0'
DESCRIPTION = 'EquilibraTor'
setup(
    name="EquilibraTor",
    version=VERSION,
    author="José D. D. Cediel-Becerra",
    author_email="jcedielbecerra@ufl.edu",
    description=DESCRIPTION,
    packages=find_packages(),
    package_data={"": ["equilibration.mdp", "equilibration_2.mdp", "ions.mdp", "minim.mdp"]},
    scripts=['bin/EquilibraTor'])
