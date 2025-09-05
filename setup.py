import setuptools
import os

# def get_requirements(path:str=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'requirements.txt')):
#     with open(path) as f:
#         requirements = f.read().splitlines()
#     return requirements

commands = list()
setuptools.setup(
    name='overdrive',
    version='0.1',    
    description='N/A',
    author='Philippa Richter',
    author_email='prichter@berkeley.edu',
    packages=['src'])
    # entry_points={'console_scripts':commands},
    # install_requires=get_requirements())
