import os
import sys
import tarfile
import zipfile
import subprocess
import logging
import platform
from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install

install_requires = [
    'six', 'setuptools >= 0.7', 'simplesam >= 0.0.3', 'tqdm', 'fastqp >= 0.2',
    'pandas', 'strandex', 'gffutils', 'pyfaidx', 'drmaa', 'twobitreader', 'requests'
]


class PostDevelopCommand(develop):
    def run(self):
        develop.run(self)
        install_binary_dependencies()


class PostInstallCommand(install):
    def run(self):
        install.run(self)
        install_binary_dependencies()


def install_binary_dependencies():
    sys.path.pop(0)  # ignore local search path
    import pisces
    import glob
    import requests
    from io import BytesIO
    from urllib.request import urlopen
    from shutil import rmtree
    from subprocess import call
    
    redist = os.path.join(os.path.dirname(pisces.__file__), 'redist')
    rmtree(os.path.join(redist, "salmon"), ignore_errors=True)
    local_redist = os.path.join(os.path.dirname(__file__), 'pisces/redist')
   
    if platform.system() == "Linux":
        salmon_url = "https://anaconda.org/bioconda/salmon/1.3.0/download/linux-64/salmon-1.3.0-hf69c8f4_0.tar.bz2"
    elif platform.system() == "Darwin":
        salmon_url = "https://anaconda.org/bioconda/salmon/1.3.0/download/osx-64/salmon-1.3.0-hb70dc8d_0.tar.bz2"

    print("Installing salmon")
    with requests.get(salmon_url, allow_redirects=True) as response:
        tar_file = BytesIO(response.content)
        tar_file.seek(0)
        with tarfile.open(fileobj=tar_file) as tar:
            tar.extractall(path=os.path.join(redist, "salmon"))

setup(
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    },
    name='novartis-pisces',
    author='Matthew Shirley',
    author_email='matt_d.shirley@novartis.com',
    url='https://github.com/Novartis/pisces',
    description=
    'PISCES: pipeline for rapid transcript quantitation, genetic fingerprinting, and quality control assessment of RNAseq libraries',
    license='Apache',
    packages=['pisces'],
    include_package_data=True,
    install_requires=install_requires,
    use_scm_version={"local_scheme": "no-local-version"},
    setup_requires=['setuptools_scm'],
    entry_points={'console_scripts': ['pisces = pisces.cli:main']},
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: Apache Software License",
        "Environment :: Console", "Intended Audience :: Science/Research",
        "Natural Language :: English", "Operating System :: Unix",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ])
