from setuptools import setup

setup(name='scipuzzle',
      version='1.1',
      description='Recreates a macrocomplex given different PDB \
      files containing interacting protein pairs.',
      keywords='macrocomplex bioinformatics structural pdb',
      url='http://github.com/gabcg/scipuzzle',
      author='Luisa Santus, Aina Rill Hinarejos, Gabriel Carbonell Gamón',
      author_email='scimatch@gmail.com',
      packages=['scipuzzle'],
      install_requires=['biopython', 'gooey'],
      include_package_data=True)
