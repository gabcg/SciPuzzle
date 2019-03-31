from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='scipuzzle',
      version='0.1',
      description='Recreates a macrocomplex given different PDB \
      files containing interacting protein pairs.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      keywords='macrocomplex bioinformatics structural pdb',
      url='http://github.com/gabcg/scipuzzle',
      author='Luisa Santus, Aina Rill Hinarejos, Gabriel Carbonell Gamón',
      author_email='scimatch@gmail.com',
      packages=['scipuzzle'],
      install_requires=['biopython', 'gooey'],
      include_package_data=True)
