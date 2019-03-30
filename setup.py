from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='scipuzzle',
      version='1.0',
      description='Recreates a macrocomplex given different PDB \
      files containing interacting protein pairs.',
      long_description=readme(),
      keywords='macrocomplex bioinformatics structural pdb',
      url='http://github.com/gabcg/scipuzzle',
      author='Luisa Santus, Aina Rill Hinarejos, Gabriel Carbonell Gamón',
      author_email='scimatch@gmail.com',
      packages=['scipuzzle'],
      install_requires=['biopython', 'gooey'],
      include_package_data=True)
