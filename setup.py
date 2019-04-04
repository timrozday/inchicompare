from setuptools import setup

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='InChi Compare',
    url='https://github.com/timrozday/inchi_compare',
    author='Tim Rozday',
    author_email='timrozday@ebi.ac.uk',
    # Needed to actually package something
    packages=['inchi_compare'],
    # Needed for dependencies
    install_requires=[],
    # *strongly* suggested for sharing
    version='0.1',
    # The license can be anything you like
    license='Do what you like with it (just nothing evil)',
    description='A few functions for splitting, normalising and comparing InChi strings.',
    # We will also need a readme eventually (there will be a warning)
    long_description='A few functions for splitting, normalising and comparing InChi strings. Detects exact matches by InChi layer which is useful for detecting ions and stereoisomers. Can split a mixture InChi into single molecule InChi keys.',
    # long_description=open('README.txt').read(),
)