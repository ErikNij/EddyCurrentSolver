from setuptools import setup, find_packages

setup(
    name='foamTools',  # The name of your package
    version='1.0.0',   # Your current version
    description='Tools for working with OpenFOAM', 
    author='Europe', 
    author_email='your.email@provider.com',
    packages=find_packages(),  # Automatically finds packages within the 'foamTools' directory
    install_requires=[  
        # List any external dependencies here (i.e., from PyPI)
    ],
)