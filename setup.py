from setuptools import setup, find_packages

setup(
    name='lomap2', 
    version='1.0.0', 
    description='Lead Optimization Mapper 2',
      
    url='https://github.com/nividic/Lomap',
    
    author='Gaetano Calabro and David Mobley',
    author_email='gcalabro@uci.edu',
    license='LGPL',
    
    packages=find_packages(),  
      
    include_package_data=True,
      
    entry_points={
        'console_scripts': [
            'lomap=lomap.dbmol:startup',
        ],
    },
)


